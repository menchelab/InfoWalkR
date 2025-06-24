# InfoWalkR
package implementation of the MultiOme Paper

```R
library(InfoWalkR, quietly = TRUE)
```

```R
#read in data
path = '<path_to_your_networks>'
graph_list = read_graphs_from_folder(path)

graph_list
```

```R
#test gene set
path_gene_set = '<path_to_your_test_gene_set>'

# Read the TSV file
data <- read.delim(path_gene_set, sep = "\t", header = TRUE)

# Extract the gene names column as a vector
# in the example, th e
gene_set <- unlist(strsplit(data$all_genes, ";"))

gene_set
```

```R
```

```R
```

```R
```

```R
```



