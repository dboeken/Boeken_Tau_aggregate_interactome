

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib

from loguru import logger
logger.info('Import OK')

font = {'family' : 'normal',
'weight' : 'normal',
'size'   : 8 }
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'


def filter_levels(summary_df, level=None, level_min=0.0, level_max=None, proteins_min=2):
    # calculate log2 fold enrichment
    summary_df['log2_fold_enrichment'] = np.log2(summary_df['fold_enrichment'])

    # assign filtering criteria to determine what level (and optionally how many proteins must be in that level) for plotting
    filtered = summary_df.copy()

    if level:
        logger.info(f'filtering to only level {level}')
        filtered = filtered[filtered['level'] == level]
    if proteins_min:
        logger.info(f'filtering min number of proteins {proteins_min}')
        filtered = filtered[filtered['number_in_list'] >= proteins_min]
    if level_min:
        logger.info(f'filtering min level {level_min}')
        filtered = filtered[filtered['level'] >= level_min]
    if level_max:
        logger.info(f'filtering max level {level_max}')
        filtered = filtered[filtered['level'] <= level_max]
    selected = filtered.groupby(['search_type', 'column', 'family'])['level'].min().reset_index()
    selected['selected'] = 'yes'
    filtered = pd.merge(filtered, selected, how='outer', on=['search_type', 'column', 'family', 'level'])
    filtered = filtered[filtered['selected'] == 'yes']
    return filtered


def plot_enrichment(filtered, colour_dict=None, filter_top=5, output_folder=False):
    # generate bar plots

    for (category, column), df in filtered.groupby(['search_type', 'column']):
        source = df.copy()
        source = source.sort_values(['log2_fold_enrichment'], ascending=False)
        if filter_top:
            source = source.iloc[:filter_top]
        fig, ax = plt.subplots(figsize=(8, source.shape[0]*0.45))
        if colour_dict:
            sns.barplot(x="log2_fold_enrichment", y="term_label", data=source, color=colour_dict[column])
        else:
            sns.barplot(x="log2_fold_enrichment", y="term_label", data=source)
        # Add term_ids
        for position, term_id in enumerate(source['term_id']):
            plt.annotate(term_id, (2, position), color='white')
        plt.xlim(0, 7)
        plt.title(f'GO enrichment for {category}', )
        plt.ylabel(None)
        plt.xlabel("Log2 Fold Enrichment")
        sns.despine(right=True, top=True)
        if output_folder:
            plt.savefig(f'{output_folder}/{column}_{category.replace(" ", "_")}.svg')
            plt.tight_layout()
            plt.savefig(f'{output_folder}/{column}_{category.replace(" ", "_")}.png')
        plt.show()


def plot_compiled_enrichment(filtered, colour_dict=None, filter_top=5, output_folder=False):
    # generate bar plots

    for column, df in filtered.groupby(['column']):
        source = df.copy()
        source = source.sort_values(['log2_fold_enrichment'], ascending=False)
        if filter_top:
            source = source.iloc[:filter_top]
        fig, ax = plt.subplots(figsize=(8, source.shape[0]*0.45))
        if colour_dict:
            sns.barplot(x="log2_fold_enrichment", y="term_label", data=source, color=colour_dict[column])
        else:
            sns.barplot(x="log2_fold_enrichment", y="term_label", data=source)
        # Add term_ids
        for position, term_id in enumerate(source['term_id']):
            plt.annotate(term_id, (2, position), color='white')
        plt.xlim(0, 7)
        plt.title('GO enrichment' )
        plt.ylabel(None)
        plt.xlabel("Log2 Fold Enrichment")
        sns.despine(right=True, top=True)
        if output_folder:
            plt.savefig(f'{output_folder}/{column}.svg')
            plt.tight_layout()
            plt.savefig(f'{output_folder}/{column}.png')
        plt.show()

if __name__ == "__main__":

    # ------------------Complete per-treatment enrichment analysis------------------
    input_path = 'results/go_enrichment/enrichment_complete.csv'
    output_folder = 'results/go_enrichment/'
    
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    # read in summary statistics
    protein_enrichment = pd.read_csv(f'{input_path}')
    protein_enrichment.drop([col for col in protein_enrichment.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

    filtered_df = filter_levels(
        protein_enrichment, level=None, level_min=0.0, level_max=None, proteins_min=1)
    
    filtered_df = filtered_df[filtered_df['column'].str.contains('Beads')].copy()

    # bubble chart - only those overrepresented
    palette = { 'up': '#a10000', 'down': '#063478', 'all': '#442178'}
    treatments = ['AT8', 'HT7', 'MC1', 'T181', 'TOMA']

    for treatment in treatments:
        for_plotting = filtered_df[filtered_df['log2_fold_enrichment'] > 0]
        for_plotting = for_plotting[for_plotting['column'].str.contains(treatment)].copy()
        x_order = for_plotting.sort_values(['pValue', 'search_type'])[
            'term_label'].unique().tolist()
        x_labels = {term_label: f'{term_id}' for term_id,
                    term_label in for_plotting[['term_id', 'term_label']].values}
        y_vals = [f'{treatment}_Beads_down', f'{treatment}_Beads_up', f'{treatment}_Beads_all']

        for_plotting['x_order'] = for_plotting['term_label'].map(
            dict(zip(x_order, range(len(x_order)))))
        for_plotting['y_order'] = for_plotting['column'].map(
            dict(zip(y_vals, range(len(y_vals)))))
        for_plotting['color'] = for_plotting['column'].str.split('_').str[-1]
        for_plotting['size'] = round(for_plotting['log2_fold_enrichment'] * 100, 2)

        fig, ax = plt.subplots(figsize=(0.35*len(x_order), 3))
        plt.grid(True)
        sns.scatterplot(
            data=for_plotting,
            x='x_order',
            y='y_order',
            hue='color',
            palette=palette,
            ax=ax,
            size='log2_fold_enrichment',
            sizes=(5, 500),
            # alpha=0.5
        )
        ax_t = ax.secondary_xaxis('top')
        ax_t.set_xticklabels(x_order, rotation=90)
        ax_t.set_xticks(range(0, len(x_order)))

        ax.set_xticks(range(0, len(x_order)))
        ax.set_xticklabels([x_labels[x] for x in x_order], rotation=90)
        plt.xlabel('')

        ax.set_yticks(range(0, len(y_vals)))
        ax.set_yticklabels(y_vals)
        plt.ylim(-0.250, 2.25)
        plt.ylabel('Change in reactivity')

        plt.legend(bbox_to_anchor=(1.0, 1.0))
        plt.savefig(f'{output_folder}{treatment}_bubblechart.svg')

