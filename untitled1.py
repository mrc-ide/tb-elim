# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 19:37:49 2024

@author: ISA
"""

import matplotlib.pyplot as plt

# Import necessary libraries
import pandas as pd

# Load the data from the Excel file
excel_path = 'Book2.xlsx'
# Read the excel file into a pandas dataframe
data = pd.read_excel(excel_path)

# Extracting the data for plotting
years = data['Year']
uk_born = data['UK born number of notifications']
non_uk_born = data['Non-UK born number of notifications']
# Using the provided code with modifications to have an x-axis label for every single year.

# Creating the plot with the y-axis starting from 0
plt.figure(figsize=(15, 8), dpi=300)
plt.plot(years, uk_born, 'o-', color='blue', label='UK Born')
plt.plot(years, non_uk_born, 'o-', color='red', label='Non-UK Born')

# Adding title and labels
plt.title('TB Notifications by Place of Birth')
plt.xlabel('Year')
plt.ylabel('Number of Notifications')

# Setting x-axis ticks to show every year
plt.xticks(years, years)

# Setting the y-axis to start from 0
plt.ylim(bottom=0)

# Adding legend
plt.legend()

# Show grid
plt.grid(True)

# Save the figure

# Display the plot
plt.show()
