
is a module that can simulate reads consisting of only
UMIs per each, or UMI+Genomic sequence per each. The general-purpose
design gives the module this name. To achieve this purpose, a case-study
CLI should look like below:

# Usage

=== "Python"

    ``` py
    import phylotres as pts

    gspl = pts.gmat.spsimseq_bulk(
        R_root='D:/Programming/R/R-4.3.2/',
        num_samples=2,
        num_genes=10,
        simulator='spsimseq',
        sv_fpn=to('data/spsimseq_bulk.h5'),
    )
    print(gspl)
    
    ```

=== "Command"

    ``` c++
    phylotres gmat_bulk \
    -rfpn D:/Programming/R/R-4.3.2/ \
    -nspl 2 \
    -ngene 10 \
    -gsimulator spsimseq \
    -wd ./phylotres/data/spsimseq_bulk.h5 \
    -is True \
    -vb True
            
    ```



# My Documentation



# Attributes
!!! Illustration

    === "Unordered List"
        
        ``` py title="bubble_sort.py"
        R_root='D:/Programming/R/R-4.3.2/',
        num_samples=2,
        num_genes=10,
        simulator='spsimseq',
        sv_fpn=to('data/spsimseq_bulk.h5'),
        
        ```

    === "Ordered List"

        ``` markdown
        1. Sed sagittis eleifend rutrum
        2. Donec vitae suscipit est
        3. Nulla tempor lobortis orci
        ```





# Output

``` py
12/12/2023 02:02:41 logger: =========>spsimseq is being used
SPsimSeq package version 1.12.0 
R[write to console]: Estimating featurewise correlations ...

R[write to console]: Selecting candidate DE genes ...

R[write to console]: Estimating densities ...

R[write to console]: Configuring design ...

R[write to console]: Simulating data ...

R[write to console]:  ...1 of 1

12/12/2023 02:02:48 logger: =========>spsimseq completes simulation
          Gene_1  Gene_2  Gene_3  Gene_4  ...  Gene_7  Gene_8  Gene_9  Gene_10
Sample_1     0.0     0.0     0.0     0.0  ...   322.0   425.0     7.0   1202.0
Sample_2     1.0     2.0     2.0     2.0  ...   633.0   423.0    87.0   1619.0

[2 rows x 10 columns]
```