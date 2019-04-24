"""Snakemake wrapper for PCA"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2018, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import matplotlib.pyplot as plt
import numpy as np
import os.path as op
import pandas as pd
import seaborn as sns
import sklearn.decomposition as skd

from snakemake.utils import makedirs

makedirs(op.dirname(snakemake.output[0]))
makedirs(op.dirname(snakemake.log[0]))

output_prefix = op.commonprefix([
    f[:-len("_PC_1_PC_2.png")]
    for f in snakemake.output
])

with open(snakemake.log[0], "w") as log:
    # Load dataset
    data = pd.read_csv(
        snakemake.input[0],
        sep="\t",
        header=0,
        index_col=0
    )

    log.write(str(data.head()))
    log.write("\n")

    # Reduce dataset to integers and floats
    data = data[list(data.select_dtypes(include=[np.number]).columns.values)]

    # Perform PCA
    # Computing components
    nbc = len(data.columns.tolist())
    skpca = skd.PCA(n_components=nbc)

    # Fitting results
    sktransform = skpca.fit_transform(data.T)
    skvar = skpca.explained_variance_ratio_
    results = pd.DataFrame(
        sktransform,
        columns=[f"PC_{i}" for i in range(1, nbc+1, 1)],
        index=data.columns.tolist()
    )
    log.write("PCA performed")
    log.write("\n")

    # Preparing plot
    sns.set(style="darkgrid")

    # Plotting results
    for k, l in zip([1, 3], [2, 4]):
        # Handle the lack of condition
        try:
            results["Conditions"] = [
                snakemake.params.conditions[i] for i in results.index
            ]
        except ValueError:
            log.write("Number of conditions: ")
            log.write(str(len(snakemake.params.conditions)))
            log.write("\nNumber of samples: ")
            log.write(str(len(results.index.tolist())))
            log.write("\n")
            raise

        # Building facetgrid
        g = sns.FacetGrid(
            results,
            hue="Conditions",
            height=13
        )

        name1 = f"PC_{k}"
        name2 = f"PC_{l}"

        # Plotting samples on the facetgrid
        g = g.map(plt.scatter, name1, name2)
        plt.title("%s (%.2f%%) and %s (%.2f%%)" % (
            name1,
            skvar[k - 1] * 100,
            name2,
            skvar[l - 1] * 100
        ))

        # Naming samples
        points_coordinates = zip(
            results.index.tolist(),
            results[name1],
            results[name2],

        )
        for label, x, y in points_coordinates:
            plt.annotate(
                label,
                xy=(x, y),
                xytext=(-5, -5),
                textcoords='offset points', ha='center', va='top',
            )

        # fitting legend
        try:
            frame = plt.legend(
                loc=snakemake.get("legend_position", "upper left"),
                frameon=True
            ).get_frame()
            frame.set_facecolor("white")
        except AttributeError:
            pass

        plt.savefig(
            f"{output_prefix}_{name1}_{name2}.png",
            bbox_inches='tight'
        )
        log.write(f"{output_prefix}_{name1}_{name2}.png saved")
        log.write("\n")
