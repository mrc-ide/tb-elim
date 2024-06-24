# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 13:46:57 2024

@author: ISA
"""

import pandas as pd

# Load the data from the Excel file
file_path = 'Book1.xlsx'
data = pd.read_excel(file_path)

# Display the first few rows of the dataframe to understand its structure
data.head()

import matplotlib.pyplot as plt

# Remove the 'Total' row for the plot
data_plot = data.iloc[:-1]

# Set the index to the time period column for better plotting
data_plot.set_index('Unnamed: 0', inplace=True)
# Adjust the plotting to include a small gap between each group of three bars corresponding to each year

# Calculate the positions for the groups of bars
group_positions = [x for x in range(len(data_plot))]

# Plot each column with a small gap between the bars within each group
fig, ax = plt.subplots(figsize=(10, 6))
bar_width = 0.2  # Smaller bar width to accommodate gaps within groups
gap = 0.05  # Gap between bars in a group

# Generate positions for each bar within the groups
for i, column in enumerate(data_plot.columns):
    bar_positions = [x + (bar_width + gap) * i for x in group_positions]
    ax.bar(bar_positions, data_plot[column], width=bar_width, label=column)

# Create the bar plot with a grid and demonstrate how to change the bar colors

# Define custom colors for the bars
colors = ['#DA70D6', '#5D3FD3', '#ffd7b7']

# Create the bar plot again without the grid

# Create the bar plot again with a faint horizontal and vertical grid

# Create the bar plot with a high DPI for ultra HD quality

fig, ax = plt.subplots(figsize=(10, 6), dpi=300)  # Set the DPI for the figure

for i, (column, color) in enumerate(zip(data_plot.columns, colors)):
    bar_positions = [x + (bar_width + gap) * i for x in group_positions]
    bars = ax.bar(bar_positions, data_plot[column], width=bar_width, label=column, color=color)
    # Add counts above the bars
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

# Set plot title and labels
ax.set_title('Access and Completion of Treatment Over Time')
ax.set_xlabel('Time Period')
ax.set_ylabel('Number of People')

# Set x-ticks to be in the middle of the grouped bars
ax.set_xticks([x + bar_width + gap for x in group_positions])
ax.set_xticklabels(data_plot.index)

# Add a faint grid to the background
ax.set_axisbelow(True)
ax.yaxis.grid(True, color='#EEEEEE')
ax.xaxis.grid(True, color='#EEEEEE')

# Place a legend outside of the plot area
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

# Make space for the legend
plt.tight_layout(rect=[0, 0, 0.85, 1])

# Show the plot
plt.show()