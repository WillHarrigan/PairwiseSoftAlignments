
# Pairwise Soft-Alignments of Protein Sequences 

### About This Repository

This repository hosts the implementation of pairwise comparisons utilizing the soft alignment
algorithm as presented in the paper "Large Language Models, Protein
Homology, Viruses, Alignments" by Harrigan et al. The primary purpose of this repository is to enable 
a single step process to easily evaluate the similarity of protein sequences in a fasta file.

This repository contains a single Python script, `sa_interface.py`
which implements soft-alignment functions to compute protein sequence similarity. The script computes the similarity between every pairwise combination of 
protein sequences in the input fasta file. During each pairwise comparison the script generates ESM-2 (33 dimesion, 650M parameter) embeddings for the protein sequences and then uses soft-alignment functions to 
compute similarity. 


## Input
```sh
python ./sa_interface.py ./protein_sequence_data/test_protein_sequences.fasta ./alignment_output.tsv
```

### Fasta File

```sh
>seq_1
MVPFYDAYGAIVAWVFSGDVETNARTGVRVKVGRGGTSFRVDLSDGLLPTVGFRKTFPKSAAAEVAWYLQGTQDATFIRKYAPLWDKFVELIDIKGGLFMEDRAVEGVKAAYGYRWRSHFGRDQIRLAVEALRKDPSDRRCYVSAWDPAEDGLGALDQRNVPCPASFTFSVLNGELHSSFFIRSSDVFVGLPYDVMGHALLMDAVAHELRLRPGIMHVTLAHAHLYESHWDLTVEMMKQEPVVPALQLPGWTLSQVERAPDDYVVRYAEEAKQLTWPAYNPRPEVVE 

>seq_2
MLIPFYDAYAAVLAWVLHGPVETNARTGVRVKVGRGGTSFRVDLSEGVLPTVGFRKTFPRSAAAEVAWYLQGTQDATFIRKYAPLWDKFVEELPSRVVGVKAAYGYRWRSHFGRDQIRLAVEALRKDPSDRRCYVSAWDPAEDGLGELGQRNVPCPAAFTFSALGEELHSSIVLRSSDVFVGLPYDVMGHALLVDAVARELGLRPGVMHVTLAHAHLYESHWDMAAEMLRQEPVVPELPLPGVALSGIEADPDGYVLSVAAEAKRHEWPSYNPKPEVVE 

>seq_3
MQGFPFPAPARDFPSLYEELLWRLMRQGSEELNERTGKRVKAYPYGGASFTLDMSGHELPVVGRRRLYPATAAAETAWYLLGTQDPTFMMRHAKVVWEKFLEDNPDQDAGASASKIIKAAYGYRWRKHFGRDQLQLAMDALDRNPSDRRVFISAWDPAEDGLGAQGQLNVPCPVGFTFSILDGRLNSTYLLRSSDVFVGLPYDVMGHALLMAAVGETLNVPLGFMTFTMAHPHIYDVHYAMADEFIMQAPVKPSILLPRWTVDQIAAEPNAYVEKVKKDGNAVPWPDFAPRPEVVQ 
```

## Output

The script outputs a .tsv file at the specified path. The file contains the following columns:

- **Score**: Soft-align score
- **RelativeScore**: Relative score
- **SimScore**: Percent Similarity score (Relative score/minimum sequence length)
- **QueryID**: First Sequence ID
- **HitID**: Second Sequence ID
- **QueryLength**: Length of the query sequence
- **HitLength**: Length of the hit sequence

### Example Output (.tsv)

```plaintext
Score   RelativeScore   SimScore    QueryID HitID   QueryLength HitLength
277     269.3876953125  0.97        seq_1   seq_2   287         279
267     251.768798828125 0.88       seq_1   seq_3   287         296
270     253.3522186279297 0.91      seq_2   seq_3   279         296
```

## Required Libraries

To run the script, you'll need to install the following libraries:

```sh
pip install Bio
pip install torch
pip install fair-esm
pip install scikit-learn
```

## Usage

1. Install the required libraries.
2. Run the script with the input FASTA file and the desired output file path.

For example:

```sh
python path/to/script/sa_interface.py path/to/fasta/fasta_sequences.fasta path/to/output/output_file.tsv
```

## Links to Libraries

- [Biopython](https://biopython.org/)
- [PyTorch](https://pytorch.org/)
- [fair-esm](https://github.com/facebookresearch/esm)
- [scikit-learn](https://scikit-learn.org/)
