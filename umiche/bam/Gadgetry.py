__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

import functools
import pathlib
import shutil


class Gadgetry:

    def __init__(
            self,
    ):
        pass

    @staticmethod
    def index_sort(
            *, # only for keywords
            do_sort: bool = True,  # False
            do_index: bool = True,  # False
            threads: int = 4,
            keep_tmp: bool = False,
            csi: bool = False,
            replace: bool = False,
    ):
        """
        wrap any function returning a BAM file path (as str or pathlib.Path)
            @bam_postprocess()
            def foo(...): return "out.bam"

        Options
        -------
        do_sort
            Whether to perform coordinate sorting (True/False)
        do_index
            Whether to build index (must be executed after sorting, True/False)
        threads
            Number of threads used for pysam.sort / pysam.index operations (integer)
        keep_tmp
            Whether to retain unsorted temporary files when replacing the original BAM (True/False)
        csi
            Generate .csi index (for ultra-large reference genomes); defaults to .bai index
        replace
            Whether to overwrite the original BAM file after sorting (True/False)
        """
        def decorator(func):
            @functools.wraps(func)
            def wrapper(*args, **kwargs):
                bam_path = pathlib.Path(func(*args, **kwargs)).resolve()

                if not bam_path.exists():
                    raise FileNotFoundError(f"Function {func.__name__} returned "
                                            f"{bam_path}, but the file does not exist.")

                # /*** ---------- sort ---------- ***/
                import pysam
                if do_sort:
                    tmp_out = bam_path.with_suffix(".unsorted.bam") if replace else bam_path
                    if replace:
                        shutil.move(bam_path, tmp_out)  # move original files

                    sorted_path = (bam_path if replace
                                   else bam_path.with_suffix(".sorted.bam"))

                    # pysam.sort receive CLI-styled params
                    pysam.sort(
                        "-@", str(threads),
                        "-o", str(sorted_path),
                        str(tmp_out if replace else bam_path)
                    )

                    if replace and not keep_tmp:
                        tmp_out.unlink()

                # /*** ---------- index ---------- ***/
                final_path = bam_path if replace or not do_sort else bam_path.with_suffix(".sorted.bam")
                if do_index:
                    index_args = ["-@", str(threads)]
                    if csi:
                        index_args.insert(0, "-c")
                    pysam.index(*index_args, str(final_path))

                return str(final_path)
            return wrapper
        return decorator
