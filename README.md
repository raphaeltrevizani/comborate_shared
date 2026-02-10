ComboRATE:

Iteratively group HLAs to lower RATE p-value

Install the dependencies with
`pip install -r requirements.txt`

USAGE:
`python comborate.py -a <allele file> -r <response file> -p <IEDB prediction file> -o <output directory>`

options:
  -h, --help            show this help message and exit
  -a, --allele-file ALLELE_FILE
                        Allele file name
  -r, --response-file RESPONSE_FILE
                        Response file name
  -o, --output OUTPUT   Output dir name

Optional:
  -p, --prediction PREDICTION
                        HLA prediction file name
  -c, --rank-cutoff RANK_CUTOFF
                        Cutoff for the binding rank (default = 25)
  -k, --response-cutoff RESPONSE_CUTOFF
                        Cutoff for the response (default = 1)


Run test with
`python comborate.py -a test/allele.txt -r test/response.txt -p test/iedb_netmhciipan_el_result.csv -o test_result`

The test results will be in the `test_result` folder created by the program in the test above and contain, for each peptide, two files:
`<peptide>_pred.csv` with the IEDB predictions isolated for the peptide
`<peptide>_rate.csv` with the RATE results for the peptide

In the `_rate.csv` file, you will also find the ComboRATE groups of restrictions, with the number of the iteration in which the group was created and the alleles that compose the group. E.g.:

The 5 first rows test_result/01_rate.csv read, for the first column:
```
2-DQA1*03:01/DQB1*03:02+DRB1*04:07+DRB1*07:01
1-DQA1*03:01/DQB1*03:02+DRB1*04:07
DQA1*03:01/DQB1*03:02
DRB1*04:07
DRB1*07:01
```

Lines 3 and 4  contain the restrictions with the two best p-values `DQA1*03:01/DQB1*03:02`, `DRB1*04:07`, which were grouped in the first iteration, making the group `1-DQA1*03:01/DQB1*03:02+DRB1*04:07`, which means
```
1- (first iteration)
DQA1*03:01/DQB1*03:02 (restriction with the lowest p-value)
DRB1*04:07 (restriction with the second lowest p-value)
```

Later, ComboRATE added `DRB1*07:01` to the 1st group, making `2-DQA1*03:01/DQB1*03:02+DRB1*04:07+DRB1*07:01`, which means
```
2- (second iteration)
DQA1*03:01/DQB1*03:02 (restriction with the lowest p-value)
DRB1*04:07 (restriction with the second lowest p-value)
DRB1*07:01 (restriction with the third lowest p-value)
```

