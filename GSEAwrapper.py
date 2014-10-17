#Rscript process.R test.txt rownames.txt c2.biocarta.v2.5.symbols.gmt phenotypes.txt result.txt
import subprocess
import os
import select
import masterdirac.models.run as r_model
import masterdirac.utils.hddata_process as hdp
import logging
import boto
from boto.exception import S3ResponseError
from boto.s3.key import Key
import datadirac.data as dd
import os.path
import pandas


def _get_source_data( working_dir, run_model ):
    """
    Downloads the data from s3 to the local machine for processing
    """
    if not os.path.exists(working_dir):
        logging.info( "Creating directory [%s]" % (
            working_dir ) )

        os.makedirs(working_dir)
    sd = run_model['source_data']
    #grab filenames we are interested in
    file_list = [f for  k, f  in sd.iteritems() if k[-4:] == 'file']
    conn = boto.connect_s3()
    bucket = conn.get_bucket( sd['bucket'] )
    for key_name in file_list:
        s3_path, fname = os.path.split(key_name)
        local_path = os.path.join(working_dir, fname)
        try:
            logging.info( "Transferring s3://%s/%s to %s" % (sd['bucket'],key_name, local_path ))
            k = Key(bucket)
            k.key = key_name
            k.get_contents_to_filename(local_path)
            logging.info("Transfer complete")
        except S3ResponseError as sre:
            logging.error("bucket:[%s] file:[%s] download." % (sd['bucket'],key_name))
            logging.error(str(sre))
            raise(sre)

def k_nearest( k, center_age, samples ):
    return [s for a,s in sorted( [((center_age - age)**2, sid) for age,sid in samples] )[:k]]

def window( start, end, samples ):
    """
    """
    my_samples = sorted(samples )
    in_window = []
    for age,sid in my_samples:
        if start <= age <= end:
            in_window.append( sid )
    return in_window

def createGSEApackage( run_id, windows ):
    """
    Generate the files for EVA

    """
    clean_run_id = '.'.join(run_id.split('-'))
    if not os.path.exists( clean_run_id ):
        os.makedirs( clean_run_id )
    run_model = r_model.get_ANRun( run_id )
    sd = run_model['source_data']
    net_config = run_model['network_config']
    #download source data
    ###DEBUG
    working_dir = os.path.join( os.getcwd(), clean_run_id )
    if not os.path.exists( working_dir ):
        os.makedirs( working_dir )
    pandas_file = os.path.join( working_dir, "expression.pnd" )
    if not os.path.exists( pandas_file ):
        _get_source_data( working_dir , run_model )
        hdg =  hdp.HDDataGen( working_dir  )
        df, _ =  hdg.generate_dataframe( run_model['source_data'], run_model['network_config'] )
        df.save( pandas_file )
    sd_obj = dd.SourceData()
    sd_obj.load_dataframe( pandas_file )
    net_table = run_model['network_config']['network_table']
    net_source = run_model['network_config']['network_source']
    sd_obj.load_net_info(net_table, net_source )

    _, meta_file = os.path.split( run_model['source_data']['meta_file'] )
    mi = dd.MetaInfo( os.path.join( clean_run_id, meta_file ) )
    strain = mi.get_strains()
    if len(strain) > 1:
        logging.warning("More than one strain, only getting first")
        logging.warning("Strains %r" % strain )
    alleles = mi.get_nominal_alleles()
    if len( alleles ) > 2:
        logging.warning("More than two alleles, only using 'WT' and other")
        logging.warning("Alleles %r" % alleles )
    if 'WT' not in alleles:
        raise Exception("Wild type not in alleles. Alleles = %r" % alleles)
    second_allele = [allele for allele in alleles if allele != 'WT'][0]
    wt_samples = mi.get_sample_ids( strain=strain[0], allele='WT' )
    comp_samples = mi.get_sample_ids( strain=strain[0], allele = second_allele)
    assert len(wt_samples) > 0
    assert len( comp_samples ) > 0
    wt_s_a = sorted( [(mi.get_age( sid), sid) for sid in wt_samples] )
    comp_s_a = sorted( [(mi.get_age( sid), sid) for sid in comp_samples] )
    comparisons = {}
    gene_names_fname = "gene_names.txt"
    with open(os.path.join(working_dir , gene_names_fname), 'w') as gnf:
        gnf.write('\n'.join(['"%s"' % gn for gn in sd_obj.source_dataframe.index]))

    logging.info("Wrote %s" % gene_names_fname )
    network_fname = "net.gmt"
    with open( os.path.join(working_dir, network_fname), 'w') as nf:
        for pw in sd_obj.get_pathways():
            nf.write( '\t'.join([pw, 'na'] + sd_obj.get_genes( pw )) + '\n' )
    logging.info("Wrote %s" % network_fname )

    for start, end in windows:
        comparisons[(start, end)] = ( window( start, end, wt_s_a), window( start, end, comp_s_a))
    result = {}
    for win, v in comparisons.iteritems():
        window_pattern = "start%iend%i" % win
        wt_s, comp_s = v
        curr_df = sd_obj.get_expression( wt_s + comp_s )
        exp_table_fname = "%s.expression.txt" % (window_pattern) 
        with open( os.path.join(working_dir, exp_table_fname), 'w') as exp:
            header = 'NAME\tDESCRIPTION\t'+'\t'.join(curr_df.columns)
            exp.write( header + '\n')
            for gene in curr_df.index:
                line = gene + '\tna' + '\t'.join(map(str, curr_df.loc[[gene]].values[0]))
                exp.write( line + '\n' )

        pheno_fname = "%s.pheno.cls" % ( window_pattern)
        with open(  os.path.join(working_dir, pheno_fname), 'w') as ph:
            n_samp = len(wt_s) + len( comp_s )
            ph.write("%i 2 1\n" % n_samp)
            ph.write("# WT MUTANT\n")
            ph.write(' '.join(['0' for s in wt_s] + ['1' for s in comp_s]) + '\n')
        params = ( exp_table_fname, gene_names_fname, network_fname, pheno_fname, "%s.%s.result" %  (clean_run_id, window_pattern ))
        params = tuple([ os.path.join(clean_run_id,p) for p in params]) 
        fin, mess = GSEA( *params )
        for m in mess:
            if len(m[1].strip()) > 0:
                logging.info("%s: %s" % (m[0], m[1]))
        """
        DEBUG
        result[win] = parse_result( params[-1] )
        #DEBUG
    t = result.keys()[0]
    n = result[t].keys()[0]

    for dt in result[t][n].keys():
        save_table( result, "%s.%s.csv" % (run_id, dt), val_type=dt )

    """
    return result

def save_table( result, file_name, val_type='pvalue' ):
    temp = sorted([(b-a, a, b) for a,b in result.keys()[:]])
    column_names = ["[%i,%i]" % (a,b) for _, a, b in temp]
    grand_table = {}
    nets = None
    for key in result.keys():
        nets = result[key].keys()
        nets.sort()
        grand_table["[%i,%i]" % key] = []
        for net in nets:
            grand_table["[%i,%i]" % key].append( result[key][net][val_type] )
    df = pandas.DataFrame.from_dict( grand_table )
    df = df[ column_names ]
    df.index = nets
    df.to_csv( file_name )



def parse_result( result_file ):
    results = []
    with open( result_file, 'r' ) as rf:
        for line in rf:
            results.append(line.strip().split())
    assert len(results) == 2
    results[1] = results[1][1:]
    assert len(results[0]) == len(results[1])
    result_dict = {}
    for label, value in zip( results[0], results[1] ):
        parsed = label.split('.')
        parsed[0] = parsed[0][1:]
        parsed[-1] = parsed[-1][:-1]
        key = '.'.join(parsed[:-1])
        if key not in result_dict:
            result_dict[key] = {}
        result_dict[key][parsed[-1]] = float( value )
    return result_dict

def GSEA( exp_mat, gene_names, pathways, phenotypes, result_dir ):
    command_string = "java -cp gsea.jar xtools.gsea.Gsea -Xmx1024m -collapse false -res %s -gmx %s -cls %s -out %s" % ( exp_mat, pathways, phenotypes, result_dir )
    messages = [('wrapper', command_string)]
    sc_p = subprocess.Popen( command_string, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True )
    reads = (sc_p.stdout, sc_p.stderr)
    cont = True
    while cont:
        cont = sc_p.poll() is None
        ret = select.select(reads, [], [])
        for fd in ret[0]:
            if fd.fileno() == sc_p.stdout.fileno():
                messages.append(('stdout', sc_p.stdout.readline().strip()))
            if fd.fileno() == sc_p.stderr.fileno():
                messages.append(('stderr', sc_p.stderr.readline().strip() ))
    line = sc_p.stdout.readline().strip()
    while line != '':
        messages.append(('stdout', line))
        line = sc_p.stdout.readline().strip()
    line = sc_p.stderr.readline().strip()
    while line != '':
        messages.append(('stderr', line))
        line = sc_p.stderr.readline().strip()
    messages.append(('wrapper', 'Complete: returned[%i]' % cont))
    return (cont, messages)

def export_csv_results( directory='.' ):
    import re
    import pandas.io.parsers as p
    import numpy as np
    windows = [(i, i+5) for i in range(4,16)] + [(4,20), (4,12), (12,20)]
    results = {}
    for d in os.listdir(directory):
        if os.path.isdir(os.path.join(directory, d)):
            run_id = d
            nets = []
            results = {}
            results['upregulated'] = {}
            results['downregulated'] = {}
            for d2 in os.listdir(os.path.join(directory,d)):
                m = re.match(r'.+start(?P<start>\d+)end(?P<end>\d+)\.result', d2)
                if m:
                    start = int( m.group('start'))
                    end = int(m.group('end'))
                    new_path = os.path.join(directory, d, d2 )
                    for d3 in os.listdir(new_path):
                        if d3[:4] == 'my_a':
                            new_path2 = os.path.join( new_path, d3 )
                            for f3 in os.listdir( new_path2 ):
                                rep =  re.match(r'gsea_report_for_(?P<type>\w+)_.*\.xls', f3)
                                if rep:
                                    if rep.group('type') != 'WT':
                                        t_type = 'downregulated'
                                    else:
                                        t_type = 'upregulated'
                                    full_path = os.path.join(directory, d,d2,d3,f3)
                                    pval_table = p.read_table(full_path)
                                    pv_dict = pval_table.set_index('NAME')['NOM p-val'].to_dict()
                                    if 'int' not in results[t_type]:
                                        results[t_type]['int'] = []
                                    results[t_type]['int'].append( (start, end) )
                                    for k,v in pv_dict.iteritems():
                                        if k not in results[t_type]:
                                            nets.append(k)
                                            results[t_type][k] = { (start, end) : v } 
                                        else:
                                            results[t_type][k][(start, end)] = v
            nets = sorted(list(set(nets)))
            up_df = pandas.DataFrame( np.ones( (len(nets), len(windows))), 
                    columns=["[%i,%i]" % ( a,b) for a,b in windows ],
                    index=nets)
            up_df.index.name = 'networks'
            down_df = up_df.copy()
            for k,v in results['upregulated'].iteritems():
                if k == 'int':
                    continue
                for inter, value in v.iteritems():
                    up_df.ix[ k ]["[%i,%i]" % inter] = value
            up_df.to_csv( "%s-upregulated-nom-pval.csv" % run_id )
            for k,v in results['downregulated'].iteritems():
                if k == 'int':
                    continue
                for inter, value in v.iteritems():
                    down_df.ix[ k ]["[%i,%i]" % inter] = value
            down_df.to_csv( "%s-downregulated-nom-pval.csv" % run_id )

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG, filename="megarun-gsea.log")
    import re
    import pandas.io.parsers as p
    import numpy as np
    print "here"
    for r in r_model.get_ANRun():
        if r['run_id'] not in ['cd1-analysis2', 'fvb-analysis1']:

            print "blah"
        else:#if r['status'] == 20 and r['run_id'][:4] not in ['test', 'lab-', 'joc-']:
            run_id = r['run_id']
            windows = [(i, i+5) for i in range(4,16)] + [(4,20), (4,12), (12,20)]
            try:
                eva_res = createGSEApackage(run_id, windows)
            except:
                logging.exception("Error running GSEA")
    export_csv_results()
