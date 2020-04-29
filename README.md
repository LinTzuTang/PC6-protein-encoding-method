# PC6-protein-encoding-method
The core idea of PC6 encoding method is using physicochemical properties as word embeddings. Each amino acid character in sequence would be replaced to a vector composed by six physicochemical property values.

# 
### Physicochemical properties clustering analysis
Given that each amino acid possesses a large number of physicochemical properties, selecting physicochemical properties that could represent the amino acids is an important step of protein encoding. Here, we collected physicochemical properties of amino acids from R package ‘Peptides’. After that, we filtered out properties that contain “NA” in the dataset and obtained the remaining 115 properties. Then, we used R function to calculate the correlation between each property and applied clustering analysis through hierarchical clustering. Finally, taking the K-means approach, we determined six as the optimal number of clusters (Supplementary Figure 2). Therefore, six physicochemical properties were chosen from the six clusters as the following: hydrophobicity (H1), volume of side chains (V), polarity (P1), pH at the isoelectric point (pl), the negative of the logarithm of the dissociation constant for the -COOH group (pKa), and net charge index of side chain (NCI). Those physicochemical properties were further selected as the features in PC6 protein encoding. 

#### fasta -> dict
```python
from Protein_Encoding import PC_6
PC_6(fasta_name)
```

![image](PC_6.png)
