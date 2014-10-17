library('GSA');
library('GSReg');
args = commandArgs(trailingOnly = TRUE);
print(args);
expression_matrix = args[1];
gene_names_file = args[2];
path_file = args[3];
phenotypes_file = args[4];
result_file = args[5];

temp_path = GSA.read.gmt(path_file);
path = list();
for (i in 1:length(temp_path$genesets) )
    path[temp_path$geneset.names[i]] = temp_path$genesets[i];
end

expression_data = as.matrix(read.table( expression_matrix , header=FALSE));
gene_names = as.vector(read.table(gene_names_file, header=FALSE)[,1]);
rownames(expression_data) = gene_names;
phenotypes = as.vector(read.table(phenotypes_file, header=0)[,1]);
#write.csv(expression_data[unlist(path['Targets_of_Hltf.A_52_P532033']),], file=paste( result_file, "Targets_of_Hltf.A_52_P532033.", sep="."))

result = GSReg.GeneSets.EVA(geneexpres=expression_data, pathways=path, 
    phenotypes=as.factor( phenotypes ) );
write.table(result, result_file);





