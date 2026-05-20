"""Report helpers for Ka/Ks calculations."""

import os
import shutil

from pgkit.core.utils import log, script_path

# ============================================================
# R visualization script
# ============================================================
def generate_r_script(output_dir):
    """Copy R script for Ka/Ks visualization"""
    # Source script location
    src_script = script_path('plot_kaks.R')

    # Destination
    dst_script = os.path.join(output_dir, 'kaks_boxplot.R')

    # Copy script
    if os.path.exists(src_script):
        shutil.copy2(src_script, dst_script)
        log(f"R script saved: {dst_script}")
    else:
        log("Warning: plot_kaks.R not found in scripts directory")

    return dst_script
