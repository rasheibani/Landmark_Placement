import csv
import json

import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import pearsonr
from scipy.stats import linregress


# convert pareto_fronts.json to pareto_fronts.csv
def json_to_csv(json_file, csv_file):
    with open(json_file) as json_file:
        data = json.load(json_file)
    with open(csv_file, 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(data[0].keys())  # write keys as header
        for row in data:
            writer.writerow(row.values())

def analyse_pareto_fronts(csv_file):


    df = pd.read_csv(csv_file)
    # Calculate ratio
    df['ratio'] = df['total_weight'] / df['sum_of_weight']

    # Groupby num_selected_vertices and calculate stats
    groups = df.groupby('num_selected_vertices')['ratio'].agg(['mean', 'median', 'min'])

    # Plot
    plt.scatter(df['num_selected_vertices'], df['ratio'], alpha=0.5)

    # Trendlines
    plt.plot(groups.index, groups['min'], label='Minimum')
    plt.plot(groups.index, groups['median'], label='Median')
    plt.plot(groups.index, groups['mean'], label='Average')

    plt.title('Trend of Uncertainty Reduction Ratio with Number of Landmarks')
    plt.xlabel('Number of Landmarks')
    plt.ylabel('Uncertainty Reduction Ratio')
    plt.legend()

    plt.savefig('plot.jpg', dpi=300)

def draw_pareto_for_one_floorplan_only(csv_file, floorplan_id):
    df = pd.read_csv(csv_file)
    df = df[df['letter'] == floorplan_id]
    plt.scatter(df['num_selected_vertices'], df['total_weight']/df['sum_of_weight'], alpha=0.5)
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







if __name__ == '__main__':
    # json_to_csv('pareto_fronts.json', 'pareto_fronts.csv')
    #analyse_pareto_fronts('pareto_fronts.csv')
    #draw_pareto_for_one_floorplan_only('pareto_fronts.csv', 'CBS_Average-Regular_Approach1')
    calculate_correlation('pareto_fronts.csv')
    print(calculate_correlation('pareto_fronts.csv'))

