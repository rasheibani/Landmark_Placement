import csv
import json

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

from scipy.stats import pearsonr
from scipy.stats import linregress

import geopandas as gpd


# convert pareto_fronts.json to pareto_fronts.csv
def json_to_csv(json_file, csv_file):
    with open(json_file) as json_file:
        data = json.load(json_file)
    with open(csv_file, 'w') as csv_file:
        writer = csv.writer(csv_file)
        headers = list(data[0].keys())  # Convert dict_keys object to a list
        writer.writerow(headers + ['all_candidates'])
        for row in data:
            # find the corresponding floorplan and all_candidates from all_candidates.csv
            with open('../data_synthetic/all_candidates.csv') as csv_file2:
                reader = csv.reader(csv_file2)
                next(reader)
                for row2 in reader:
                    if row['letter'] == row2[0]:
                        row['all_candidates'] = int(row2[1]) / 2
                        break
            writer.writerow(row.values())


import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker

def analyse_pareto_fronts(csv_file, all_candidates=False, Complexity = False):
    df = pd.read_csv(csv_file)
    df['ratio'] = df['total_weight'] / df['sum_of_weight']

    # filter out if the num_selected_vertices is more than 5 but the ratio is less than 0.4
    df = df[~((df['num_selected_vertices'] > 5) & (df['ratio'] < 0.4))]

    if all_candidates:
        df['num_selected_vertices'] = df['num_selected_vertices'] / df['all_candidates'] * 100

    if Complexity:
        color_map = {'blue':0, 'green':1, 'yellow':2, 'orange':3, 'red':4}
        # plot the scatter in a way that the color of the points is based on the complexity
        # in a way that complexity of 0 to 0.2 is blue, 0.2 to 0.4 is green, 0.4 to 0.6 is yellow,
        # 0.6 to 0.8 is orange, 0.8 to 1 is red
        # create a new column for the color
        df['color'] = pd.cut(df['Normalized Graph Asymmetry'], bins=[0, 0.2, 0.4, 0.6, 0.8, 1], labels=['blue', 'green', 'yellow', 'orange', 'red'])
        # plot the scatter
        plt.scatter(df['num_selected_vertices'], df['ratio']*100, alpha=0.5, c=df['color'].map(color_map))
    else:
        plt.scatter(df['num_selected_vertices'], df['ratio']*100, alpha=0.2)

    if all_candidates:
        groups = df.groupby(df['num_selected_vertices'].astype(int))['ratio'].mean()
        plt.plot(groups.index, groups.values * 100, label='Mean', color='red')
        plt.xlabel('Percentage of Selected Landmarks from All Candidates')
    else:
        groups = df.groupby('num_selected_vertices')['ratio'].mean()
        plt.plot(groups.index, groups.values * 100, label='Mean', color='red')
        plt.xlabel('Number of Landmarks')

    plt.title('Trend of Uncertainty Reduction Ratio with Number of Landmarks')
    plt.ylabel('Uncertainty Reduction Ratio')

    # Show y_ticks as percentage
    plt.gca().yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())
    plt.gca().xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())

    # make sure every tick is shown with 1% interval
    plt.gca().yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=10))
    plt.gca().xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=1))

    plt.legend()
    plt.savefig('plot.jpg', dpi=300)





def draw_pareto_for_one_floorplan_only(csv_file, floorplan_id, all_candidates=False, Extended=False):
    df = pd.read_csv(csv_file)
    df = df[df['letter'] == floorplan_id]
    if all_candidates:
        df['num_selected_vertices'] = df['num_selected_vertices'] / df['all_candidates']
    if Extended:
        plt.scatter(df['num_selected_vertices'], df['total_weight'] / df['sum_of_weight'], alpha=0.5)
    plt.scatter(df['num_selected_vertices'], df['total_weight'] / df['sum_of_weight'], alpha=0.5)
    plt.title('Pareto Front')
    plt.xlabel('Number of Landmarks')
    plt.ylabel('Uncertainty Reduction Ratio')
    plt.savefig('pareto_front_' + floorplan_id + '.jpg', dpi=300)


# calculate correlation BETWEEN NUMBER OF LANDMARKS AND UNCERTAINTY REDUCTION RATIO
def calculate_correlation(csv_file):
    df = pd.read_csv(csv_file)
    df['ratio'] = df['total_weight'] / df['sum_of_weight']
    # calculate Correlation Coefficient (r) and p-value
    correlation = df['num_selected_vertices'].corr(df['ratio'])

    # Calculate Pearson correlation coefficient and p-value
    correlation, p_value = pearsonr(df['ratio'], df['num_selected_vertices'])

    # Number of data points
    n = len(df)

    # Degrees of freedom
    df_regression = 2 - 1
    df_residual = n - 2

    # Critical p-value for two-tailed test
    alpha = 0.05
    critical_p_value = alpha / 2

    # Calculate F-value and critical F-value for regression
    slope, intercept, r_value, p_value_regression, std_err = linregress(df['ratio'], df['num_selected_vertices'])
    y_predicted = slope * df['ratio'] + intercept
    sse = ((df['num_selected_vertices'] - y_predicted) ** 2).sum()
    sst = ((df['num_selected_vertices'] - df['num_selected_vertices'].mean()) ** 2).sum()
    msr = sst / df_regression
    mse = sse / df_residual
    f_value = msr / mse
    critical_f_value = 3.840  # F-value with alpha = 0.05 and df_regression = 1, df_residual = n - 2

    return correlation, p_value, critical_p_value, f_value, critical_f_value


def count_all_candidates_for_all():
    floorplan_Bbox = gpd.read_file("../data_synthetic/ICD_CalculatedV7.shp")
    floorplan_corner = gpd.read_file("../data_synthetic/CornersV7.shp")

    to_write = []

    for index, row in floorplan_Bbox.iterrows():
        polygon = row['geometry']
        letter = row['distinct']
        floorplan_Corners = floorplan_corner[floorplan_corner.intersects(polygon)]
        all_candidate = len(floorplan_Corners)
        to_write.append([letter, all_candidate])

    with open('../data_synthetic/all_candidates.csv', 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['letter', 'all_candidates'])
        for row in to_write:
            writer.writerow(row)


def add_all_candidate_column_to_csv(csv_file, all_candidates_file):
    df = pd.read_csv(csv_file)
    all_candidates = pd.read_csv(all_candidates_file)
    df = df.merge(all_candidates, on='letter')
    df.to_csv(csv_file, index=False)

def add_complexityMeasures_from_xlsx(csv_file, xlsx_file):
    df = pd.read_csv(csv_file)
    complexity = pd.read_excel(xlsx_file)
    df = df.merge(complexity, on='letter')
    df.to_csv(csv_file, index=False)


import pandas as pd
import matplotlib.pyplot as plt


def draw_subplots_based_on_complexity(csv_file, Complexity):
    df = pd.read_csv(csv_file)
    df['ratio'] = df['total_weight'] / df['sum_of_weight']

    # Define color map and corresponding colors
    color_map = {'0.0-0.25': '#b2182b', '0.25-0.50': '#fddbc7', '0.50-0.75': '#67a9cf', '0.75-1.0': '#2166ac'}

    # Create 'color' column based on complexity bins
    df['color'] = pd.cut(df[Complexity], bins=[0, 0.25, 0.50, 0.75, 1.0],
                         labels=['0.0-0.25', '0.25-0.50', '0.50-0.75', '0.75-1.0'])

    fig, axs = plt.subplots(2, 2, figsize=(10, 8), sharey=True)  # 2 rows, 2 columns

    for (key, color), ax in zip(color_map.items(), axs.ravel()):
        group = df[df['color'] == key]
        # filter out if the num_selected_vertices is more than 5 but the ratio is less than 0.4
        group = group[~((group['num_selected_vertices'] > 5) & (group['ratio'] < 0.4))]

        ax.scatter(group['num_selected_vertices'], group['ratio'] * 100, label=key, alpha=0.7, color=color)
        # make sure that the x-axis starts from 0 end at 15
        ax.set_xlim(0, 15)


        #draw A dashed vertical line and show the value on x-axis point fot X: mean of all values of num_selected_vertices
        # make sure the x-axis value is shown where the vertical line is drawn
        ax.axvline(group['num_selected_vertices'].mean(), linestyle='--', color='black', alpha=0.7)
        ax.text(group['num_selected_vertices'].mean() + 0.1, 0.2*100, f'{group["num_selected_vertices"].mean():.2f}', verticalalignment='center')



        ax.set_title(f'Complexity: {key}')
        ax.set_xlabel('Number of Landmarks')
        ax.set_ylabel('Uncertainty Reduction Ratio')

        # Calculate statistical measure (mean, median, etc.)
        statistical_measure = group['num_selected_vertices'].mean()  # Change this to the desired statistical measure

        # Add statistical measure as text to the plot
        # ax.text(0.3751, 0.15, f'mean of required landmarks: {statistical_measure:.2f}', transform=ax.transAxes, verticalalignment='top',
        #         bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.5))

    plt.tight_layout()  # Adjust the layout to prevent overlap
    plt.savefig('plotsubplots_with_statistical_measure.jpg', dpi=300)

def draw_heatmap(csv_file):
    df = pd.read_csv(csv_file)
    df['ratio'] = df['total_weight'] / df['sum_of_weight']
    # filter out if the num_selected_vertices is more than 5 but the ratio is less than 0.4
    df = df[~((df['num_selected_vertices'] > 5) & (df['ratio'] < 0.4))]

    # Draw a 2D density plot
    plt.hexbin(df['num_selected_vertices']/df['all_candidates']*100, df['ratio']*100, gridsize=13, cmap='Blues', bins='log')
    plt.colorbar(label='log10(count)')
    plt.title('Heatmap of Uncertainty Reduction Ratio')
    plt.xlabel('Percentage of Landmarks Selected from all Candidates')
    plt.ylabel('Uncertainty Reduction Ratio (%)')
    plt.savefig('heatmap.jpg', dpi=300)





def draw_aggregation_of_all_results2(csv_files):
    fig, axs = plt.subplots(4, 4, figsize=(16, 16), sharey=True)

    color_map = {'0.0-0.25': '#b2182b', '0.25-0.50': '#fddbc7', '0.50-0.75': '#67a9cf', '0.75-1.0': '#2166ac'}

    for i, csv_file in enumerate(csv_files):
        df = pd.read_csv(csv_file)

        df['ratio'] = df['total_weight'] / df['sum_of_weight']

        df = df[~((df['num_selected_vertices'] > 5) & (df['ratio'] < 0.4))]

        df['color'] = pd.cut(df['Normalized Graph Asymmetry'], bins=[0, 0.25, 0.50, 0.75, 1.0],
                             labels=['0.0-0.25', '0.25-0.50', '0.50-0.75', '0.75-1.0'])

        counter = 0

        for color_range in ['0.0-0.25', '0.25-0.50', '0.50-0.75', '0.75-1.0']:
            ax = axs[counter, i]
            # set the font size of all texts in the plot
            plt.rcParams.update({'font.size': 14})
            group = df[df['color'] == color_range]
            ax.scatter(group['num_selected_vertices'], group['ratio'] * 100, label=color_range, alpha=1, color=color_map[color_range])
            ax.axvline(group['num_selected_vertices'].median(), linestyle='--', color='black', alpha=0.7)
            ax.text(group['num_selected_vertices'].median(), 0.1*100, f'{group["num_selected_vertices"].mean():.2f}', verticalalignment='top', horizontalalignment='left')
            ax.set_xlim(0, 15)
            # do not show top and right axis
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)


            if counter == 3:
                ax.set_xlabel('Number of Landmarks', fontsize=14)
            else:
                ax.set_xticklabels([])

            if i == 0:
                ax.set_ylabel(f'(Complexity: {color_range})\n'+'Uncertainty Reduction Ratio', fontsize=14)

            if counter==0:
                if i==0:
                    ax.set_title('4 Sectors')
                elif i==1:
                    ax.set_title('6 Sectors')
                elif i==2:
                    ax.set_title('8 Sectors')
                else:
                    ax.set_title('Klippel')
            counter += 1

    axs[3, 0].set_xlabel('Number of Landmarks')
    axs[3, 1].set_xlabel('Number of Landmarks')



    # save the figure in with a ration of 2 to 1 and 300 dpi
    plt.tight_layout()
    plt.savefig('aggreagte2.jpg', dpi=300)




if __name__ == '__main__':

    # json_to_csv('pareto_fronts.json', 'pareto_fronts.csv')
    # count_all_candidates_for_all()
    # add_all_candidate_column_to_csv('pareto_fronts.csv', 'all_candidates.csv')
    # add_complexityMeasures_from_xlsx('pareto_fronts.csv', 'Complexity.xlsx')
    # draw_pareto_for_one_floorplan_only('pareto_fronts.csv', 'CBS_Average-Regular_Approach1', all_candidates=False)
    # draw_subplots_based_on_complexity('pareto_fronts.csv', 'Normalized Graph Asymmetry')
    # analyse_pareto_fronts('pareto_fronts.csv', all_candidates=True, Complexity=False)
    # draw_heatmap('pareto_fronts.csv')
    # print(calculate_correlation('pareto_fronts.csv'))
    draw_aggregation_of_all_results2(['V2results-4sectors/pareto_fronts.csv', 'V2results-6sectors/pareto_fronts.csv', 'V2results-8sectors/pareto_fronts.csv', 'V2results-Klippel/pareto_fronts.csv'])