The deduplicated UMI counts are saved uniformly in `{method}_dedup.txt`.

The column name means the round of permutation, e.g. `pn0` means the 1st round of permutation with respect to a certain sequencing condition, e.g. sequencing errors from 0.000001 to 0.1.

``` shell
pn0	pn1	pn2	pn3	pn4	pn5	pn6	pn7	pn8	pn9
50	50	50	50	50	50	50	50	50	49
50	50	50	50	50	50	50	50	50	50
50	50	50	50	50	50	50	50	50	50
50	50	50	50	50	50	50	50	50	50
51	51	51	51	51	51	51	51	51	51
51	51	51	51	51	51	51	51	51	51
51	51	51	51	51	51	51	52	51	51
51	51	51	51	51	51	52	51	51	51
50	50	50	50	50	50	50	50	50	50
51	51	51	51	51	51	51	51	51	51
50	50	51	50	51	51	50	51	51	50
51	51	51	51	51	51	51	51	51	51
51	51	51	51	51	51	51	51	51	51

```


In `set_cover_multi_len_split_by_mv_dedup.txt`, each cell means in each condition (e.g. at sequencing error rate 0.001) how many times are needed for set covering UMIs, with each element separated by `;` meaning how many steps are needed to merge UMIs each time.

``` shell
pn0	pn1	pn2	pn3	pn4	pn5	pn6	pn7	pn8	pn9
1;1;1;1	1;1;1;1	1;1;1	1;1;1	1;1;1	1;1	1;1;1	1;1;1	1;1;1;1	1;1;1
1;1;1;1;1;1	1;1;1;1;1;1	1;1;1;1;1	2;1;1;1	2;1;1;1	1;1;1;1	1;1;1;1;1	1;1;1;1;1	1;1;1;1;1;1	1;1;1;1;1
2;1;1;1;1;1;1;1;1;1;1;1	2;1;1;1;1;1;1;1;1;1;1;1	2;1;1;1;1;1;1;1;1;1;1	2;2;1;1;1;1;1;1;1;1	2;2;1;1;1;1;1;1;1;1	2;1;1;1;1;1;1;1;1;1	2;1;1;1;1;1;1;1;1;1;1	
...
138;130;126;126;124;120;119;117;115;114;112;112;111;110;110;110;110;109;109;107;107;106;106;105;103;102;100;100;99;97;91;91;90;89;89;86;85;75;73;70;70;63;63;57;56;45;44;41;35;27;5;5;5;5;5;4;4;4;4;4;4;4;4;4;4;4;4;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1
```

Below shows how many unique UMIs are or are not solved by set cover. Their combination means the total count of unique UMIs after set cover deduplication.

=== ":material-file-table-outline:`set_cover_solved_split_by_mv_dedup.txt`"

    ``` shell
    pn0	pn1	pn2	pn3	pn4	pn5	pn6	pn7	pn8	pn9
    4.0	4.0	3.0	3.0	3.0	2.0	3.0	3.0	4.0	3.0
    6.0	6.0	5.0	4.0	4.0	4.0	5.0	5.0	6.0	5.0
    12.0	12.0	11.0	10.0	10.0	10.0	11.0	10.0	12.0	11.0
    16.0	16.0	14.0	14.0	13.0	14.0	15.0	14.0	16.0	15.0
    20.0	19.0	18.0	18.0	17.0	18.0	19.0	18.0	20.0	19.0
    32.0	31.0	31.0	31.0	30.0	31.0	31.0	31.0	33.0	32.0
    40.0	39.0	40.0	40.0	39.0	39.0	39.0	39.0	39.0	40.0
    48.0	48.0	48.0	48.0	48.0	48.0	48.0	48.0	48.0	48.0
    49.0	49.0	49.0	49.0	49.0	49.0	49.0	49.0	49.0	49.0
    50.0	50.0	50.0	50.0	50.0	50.0	50.0	50.0	50.0	50.0
    50.0	50.0	50.0	50.0	50.0	50.0	50.0	50.0	50.0	50.0
    50.0	50.0	50.0	50.0	50.0	50.0	50.0	50.0	50.0	50.0
    50.0	50.0	50.0	50.0	50.0	50.0	50.0	50.0	50.0	50.0
    52.0	51.0	52.0	51.0	51.0	51.0	52.0	53.0	50.0	51.0
    64.0	63.0	60.0	64.0	60.0	62.0	57.0	61.0	55.0	62.0
    91.0	95.0	97.0	97.0	88.0	104.0	88.0	76.0	91.0	87.0
    155.0	154.0	144.0	146.0	161.0	143.0	168.0	156.0	154.0	149.0
    519.0	544.0	571.0	558.0	522.0	551.0	539.0	533.0	555.0	557.0
    ```

=== ":material-file-table-outline:`set_cover_not_solved_split_by_mv_dedup.txt`"

    ``` shell
    pn0	pn1	pn2	pn3	pn4	pn5	pn6	pn7	pn8	pn9
    46.0	46.0	47.0	47.0	47.0	48.0	47.0	47.0	46.0	47.0
    44.0	44.0	45.0	46.0	46.0	46.0	45.0	45.0	44.0	45.0
    38.0	38.0	39.0	40.0	40.0	40.0	39.0	40.0	38.0	39.0
    34.0	34.0	36.0	36.0	37.0	36.0	35.0	36.0	34.0	35.0
    30.0	31.0	32.0	32.0	33.0	32.0	31.0	32.0	30.0	31.0
    18.0	19.0	19.0	19.0	20.0	19.0	19.0	19.0	17.0	18.0
    10.0	11.0	10.0	10.0	11.0	11.0	11.0	11.0	11.0	10.0
    2.0	2.0	2.0	2.0	2.0	2.0	2.0	2.0	2.0	2.0
    1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0
    0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
    3.0	3.0	1.0	1.0	1.0	2.0	1.0	0.0	1.0	2.0
    4.0	6.0	4.0	4.0	6.0	3.0	2.0	2.0	2.0	4.0
    8.0	7.0	10.0	8.0	5.0	9.0	9.0	7.0	7.0	8.0
    49.0	51.0	56.0	50.0	55.0	44.0	51.0	59.0	56.0	53.0
    153.0	149.0	175.0	172.0	169.0	163.0	161.0	165.0	165.0	160.0
    328.0	331.0	328.0	313.0	336.0	334.0	343.0	339.0	334.0	358.0
    467.0	455.0	463.0	472.0	482.0	483.0	451.0	460.0	516.0	501.0
    1002.0	975.0	918.0	933.0	960.0	933.0	995.0	969.0	964.0	919.0

    ```

# Read UMI trajectory files

The UMI trajectories across PCR amplification cycles are represented to be in a Python multiple dimensional dictionary. From outer to inner, the keys stand for:

**Permutation number** :material-arrow-right-thin: **sequencing condition** :material-arrow-ri  ght-thin: connected component :material-arrow-right-thin: **UMI that are merging other UMIs** :material-arrow-right-thin: tuple (1. **origin** (UMI identity before PCR amplification); 2. **merged UMIs that have the same origin**; 3. **merged UMIs that originate differently**.

``` shell
{"0": {"1e-05": {"cc0": {}, "cc1": {}, "cc2": {}, "cc3": {}, "cc4": {}, "cc5": {}, "cc6": {}, "cc7": {}, "cc8": {}, "cc9": {}, "cc10": {}, "cc11": {}, "cc12": {}, "cc13": {}, "cc14": {}, "cc15": {}, "cc16": {}, "cc17": {}, "cc18": {}, "cc19": {}, "cc20": {}, "cc21": {}, "cc22": {}, "cc23": {}, "cc24": {}, "cc25": {}, "cc26": {}, "cc27": {}, "cc28": {}, "cc29": {}, "cc30": {}, "cc31": {}, "cc32": {}, "cc33": {}, "cc34": {}, "cc35": {}, "cc36": {}, "cc37": {}, "cc38": {}, "cc39": {}, "cc40": {}, "cc41": {}, "cc42": {}, "cc43": {}, "cc44": {}, "cc45": {}, "cc46": {}, "cc47": {}, "cc48": {}, "cc49": {}}, "2.5e-05": {"cc0": {}, "cc1": {}, "cc2": {}, "cc3": {}, "cc4": {}, "cc5": {}, "cc6": {}, "cc7": {}, "cc8": {}, "cc9": {}, "cc10": {}, "cc11": {}, "cc12": {}, "cc13": {}, "cc14": {}, "cc15": {}, "cc16": {}, "cc17": {}, "cc18": {}, "cc19": {}, "cc20": {}, "cc21": {}, "cc22": {}, "cc23": {}, "cc24": {}, "cc25": {}, "cc26": {"26": {"ori": 41, "same": [50, 26], "diff": []}}, "cc27": {}, "cc28": {}, "cc29": {}, "cc30": {}, "cc31": {}, "cc32": {}, "cc33": {}, "cc34": {}, "cc35": {}, "cc36": {}, "cc37": {}, "cc38": {}, "cc39": {}, "cc40": {}, "cc41": {}, "cc42": {}, "cc43": {}, "cc44": {}, "cc45": {}, "cc46": {}, "cc47": {}, "cc48": {}, "cc49": {}}, "5e-05": {"cc0": {}, "cc1": {}, "cc2": {}, "cc3": {}, "cc4": {}, "cc5": {}, "cc6": {}, "cc7": {}, "cc8": {}, "cc9": {}, "cc10": {}, "cc11": {}, "cc12": {}, "cc13": {}, "cc14": {}, "cc15": {}, "cc16": {}, "cc17": {}, "cc18": {}, "cc19": {}, "cc20": {}, "cc21": {}, "cc22": {}, "cc23": {}, "cc24": {}, "cc25": {"25": {"ori": 25, "same": [51, 25], "diff": []}},
```