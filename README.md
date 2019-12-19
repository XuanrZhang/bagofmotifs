# bagofmotifs

## Description 
Program searches for orthologous gene enhancers between query and target sequences using transcription factor binding site motif signatures. 
<Br>
<Br>
Aligns divergent enhancers based on their collection of binding motifs rather than the standard approach of using sequences alone.

## Usage

To run the program from the command line, enter:

```
./main.R <query file> <target file>
```

- Relative path to query file (lasagna motif file)
- Relative path to target file (lasagna motif file)

## Command line arguments

The first two arguments must be the query and target motif files.

The following optional flags can be used to set parameters in the program 

Flag | Short| Type | Constraints| Description | Default
---- | --- | --- | --- | --- | ---
-pval | -p | double |0 <= pval <= 1| Motif discovery threshold | 0.001 
-window | -w |integer | window >= 1| Window size for sliding window function | 900
-seed | -s | integer |any integer| Using the same seed will ensure the same sample is taken when calculating emperical null distribution | NULL

## Input file format

## Output format

By default output is directed to STDOUT and can be redirected into desired file. e.g.:
```
./main.R <query file> <target file> > output.txt
```






