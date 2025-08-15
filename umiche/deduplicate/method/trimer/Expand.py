class Expand:

    def __init__(
            self,
    ):
        # The IUPAC legal characters (including common ambiguous bases) are:
        # A, C, G, T, U, R, Y, S, W, K, M, B, D, H, V, N.
        self.IUPAC = set("ACGTNRYKMSWBDHV")

    def homotrimer_umi(
            self,
            umi: str,
            validate: bool = True,
    ) -> str | None:
        """
        """
        if umi is None:
            return None
        s = str(umi).strip().upper()
        if not s:
            return None

        if validate:
            if len(s) != 8:
                raise ValueError(f"Expected length-8 UMI, got {len(s)}: {umi!r}")
            bad = set(s) - self.IUPAC
            if bad:
                raise ValueError(f"UMI contains invalid base(s): {''.join(sorted(bad))} in {umi!r}")

        return "".join(base * 3 for base in s)


if __name__ == "__main__":
    p = Expand()
    ht_umi = p.homotrimer_umi("ACGTACGT")
    print(ht_umi)

    import pandas as pd
    df = pd.read_csv(
        '/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/umicountr/df_res.txt',
        sep='\t',
        header=0,
    )
    df["htUMI"] = df["UX"].apply(p.homotrimer_umi)
    print(df)