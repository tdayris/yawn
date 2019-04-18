#!/usr/bin/python3
# conding: utf-8

"""
Plot a PCA on selected samples

Unlincense terms of use:
This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <http://unlicense.org/>
"""

import begin
from collections import OrderedDict
import colorsys
import itertools
import logging
import matplotlib.pyplot as plt
import multiqc.plots as mup
import numpy as np
import pandas as pd
import seaborn as sns
import sklearn.decomposition as skd
import traceback
import yaml


@begin.convert(paths=str)
def readSalmon(*paths: "Multiple paths to Salmon files"):
    """Return a DataFrame from salmon files"""
    frames = None
    for path in paths:
        logging.debug("Reading Salmon file: %s" % path)
        data = pd.read_csv(path, sep="\t", header=0, index_col=0)
        data = data[["TPM"]]
        data.columns = [path]
        try:
            frames = pd.merge(frames, data, left_index=True, right_index=True)
        except ValueError:
            frames = data

    return frames


@begin.start
@begin.logging
@begin.tracebacks
@begin.convert(salmon_paths=str, title=str, output=str, legend_position=str,
               conditions=str, show=begin.utils.tobool, drop_sample=str)
def main(normalized_table: "Path to TPM file",
         *conditions: "Space separated list of conditions per sample",
         drop_sample: "Comma separated list of samples to be dropped out" = "",
         legend_position: "Position of the legend" = "upper center",
         output: "PCA output prefix" = "PCA",
         samples_names: "Write sample name aside each point" = True,
         show: "Should the graph be displayed instead of saved ?" = False,
         multiqc_yaml_config: "Output copy of results "
                              "in Yaml for MultiQC" = False) -> None:
    """Plot a PCA on selected samples"""
    # Load datasets
    data = pd.read_csv(normalized_table, sep="\t", header=0, index_col=0)
    logging.debug(data.head())
    data = data[list(data.select_dtypes(include=[np.number]).columns.values)]
    logging.debug(data.head())

    condition_dict = {
        k: v for k, v in zip(data.columns, conditions)
    }

    if "," in drop_sample:
        data = data.T[~data.columns.isin(drop_sample.split(","))].T

    # Perform pca
    logging.debug("Computing components")
    nbc = len(data.columns.tolist())
    skpca = skd.PCA(n_components=nbc)

    # Prepare plots
    logging.debug("Fitting results")
    sktransform = skpca.fit_transform(data.T)
    logging.debug(data.T.head())
    skvar = skpca.explained_variance_ratio_
    results = pd.DataFrame(
        sktransform,
        columns=["PC_%i" % i for i in range(1, nbc+1, 1)],
        index=data.columns.tolist()
    )

    if multiqc_yaml_config is False:
        logging.debug("Building graph")
        sns.set(style="darkgrid")
        # Handle the lack of condition
        for k, l in zip([1, 3], [2, 4]):
            try:
                results["Conditions"] = [
                    condition_dict[i] for i in results.index
                ]
            except ValueError:
                print(len(conditions), len(results.index.tolist()))
                print(conditions, results.index.tolist())
                raise
            g = sns.FacetGrid(
                results,
                hue="Conditions",
                size=13
            )

            name1 = "PC_%i" % k
            name2 = "PC_%i" % l

            g = g.map(plt.scatter, name1, name2)
            plt.title("%s (%.2f%%) and %s (%.2f%%)" % (
                name1,
                skvar[k - 1] * 100,
                name2,
                skvar[l - 1] * 100
            ))

            if samples_names:
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

            try:
                frame = (plt.legend(loc=legend_position, frameon=True)
                            .get_frame())
                frame.set_facecolor("white")
            except AttributeError:
                pass

            logging.debug("Plotting to %s" % (
                "standard output" if show else output)
            )
            if show:
                plt.tight_layout()
                plt.show()
            else:
                plt.savefig(
                    "%s_%s_%s.png" % (output, name1, name2),
                    bbox_inches='tight'
                )
    else:
        try:
            zipped = zip(
                set(conditions),
                sns.color_palette("Blues", len(set(conditions)))
            )
            colors = {k: v for k, v in zipped}
            results["color"] = [
                '#%02x%02x%02x' % (j, k, l) for j, k, l in [
                    list(map(
                        lambda x: int(round(x*255)),
                        colorsys.hls_to_rgb(*i))
                    )
                    for i in [colors[v] for v in conditions]]
            ]
        except ValueError:
            print(conditions, results.index.tolist())
            raise

        for k, l in zip([1, 3], [2, 4]):
            logging.debug("Building MultiQC output for PC%s and PC%s" % (k, l))
            cols_rename = {"PC_%s" % k: "x", "PC_%s" % l: "y"}
            sub = results[["PC_%s" % k, "PC_%s" % l, "color"]]
            data = sub.rename(columns=cols_rename).T.to_dict()

            yaml_dict = {
                "id": "my_pca_%s%s_section" % (k, l),
                "section_name": "PCA Analysis (Axies %s and %s)" % (k, l),
                "description": ("This plot shows components "
                                "from a Principal Component Analysis"),
                "plot_type": "scatter",
                "pconfig": {
                    "id": "pca_12_scatter_plot",
                    "xlab": "PC%s" % k,
                    "ylab": "PC%s" % l
                },
                "data": data
            }

            with open("%s_%s%s_mqc.yaml" % (output, k, l), "w") as yout:
                content = yaml.dump(yaml_dict, default_flow_style=False)
                yout.write(content)

        logging.debug("Building histogramm of PCA loadings")
        yaml_dict = {
            "id": "pca_loadings",
            "section_name": "PCA Loadings",
            "description": "The variance exmplained by each component",
            "plot_type": "bargraph",
            "pconfig": {
                "id": "pca_loadings",
                "xlab": "Components",
                "ylab": "Percentage of variance explained"
            },
            "data": {
                k: {"Explained": v, "Lost": 1 - v}
                for k, v in zip(
                    ["PC_%i" % i for i in range(1, len(skvar)+1)],
                    skvar
                )
            }
        }
        logging.debug(yaml_dict)

        with open("%s_pca_loadings_mqc.yaml" % (output), "w") as yout:
            content = yaml.dump(yaml_dict, default_flow_style=False)
            yout.write(content)
