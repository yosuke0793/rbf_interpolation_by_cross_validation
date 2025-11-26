import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 14
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5
plt.rcParams['axes.formatter.limits'] = (-2, 3)
plt.rcParams['axes.formatter.use_mathtext'] = True

def read_data(file_path1, file_path2):
    """
    Reads true and predicted values from a given file.

    Parameters:
    file_path (str): The path to the data file.

    Returns:
    tuple: A tuple containing two numpy arrays - true values and predicted values.
    """
    data = np.loadtxt(file_path1)
    inputs = data[:, 0:6]
    true_values = data[:, 6:9]
    predicted_values = data[:, 9:12]

    data = np.loadtxt(file_path2)
    errors = data[:, 6:9]
    return inputs, true_values, predicted_values, errors

def Rcoeff(true_values, predicted_values):
    """
    Calculates the correlation coefficient (R) between true and predicted values.

    Parameters:
    true_values (array-like): The true values.
    predicted_values (array-like): The predicted values from the model.

    Returns:
    float: The correlation coefficient R.
    """
    true_mean = np.mean(true_values, axis=0)
    predicted_mean = np.mean(predicted_values, axis=0)

    numerator = np.sum((true_values - true_mean) * (predicted_values - predicted_mean), axis=0)
    denominator = np.sqrt(np.sum((true_values - true_mean) ** 2, axis=0) * np.sum((predicted_values - predicted_mean) ** 2, axis=0))

    R = numerator / denominator
    return R

def visualize_results(inputs, true_values, predicted_values, errors, CRs):
    """
    Visualizes the true vs predicted values using a scatter plot and a line plot.

    Parameters:
    inputs (array-like): The input values.
    true_values (array-like): The true values.
    predicted_values (array-like): The predicted values from the model.
    errors (array-like): The errors between true and predicted values.
    """
    fig = plt.figure(figsize=(18, 12))
    fig.subplots_adjust(wspace=0.3, hspace=0.4)
    titleposition = -0.2
    titles = ['(a) True vs Predicted, Stress 11', 
              '(b) True vs Predicted, Stress 22',
              '(c) True vs Predicted, Stress 12',
              '(d) Errors vs True, Stress 11',
              '(e) Errors vs True, Stress 22',
              '(f) Errors vs True, Stress 12']
    axes = []
    for i in range(6):
        ax = fig.add_subplot(2, 3, i + 1)
        axes.append(ax)
    
    # Scatter plot of true vs predicted values
    for i in range(3):
        axes[i].scatter(true_values[:,i], predicted_values[:,i], color='blue', label='Predicted vs True', alpha=0.6)
        axes[i+3].scatter(true_values[:,i], errors[:,i], color='blue', label='Errors', alpha=0.6)
        axes[i].text(0.05, 0.9, f'R = {CRs[i]:.5e}', transform=axes[i].transAxes, fontsize=12, verticalalignment='top')
    
        # Line plot for perfect prediction
        min_val = min(min(true_values[:,i]), min(predicted_values[:,i]))
        max_val = max(max(true_values[:,i]), max(predicted_values[:,i]))
        axes[i].plot([min_val, max_val], [min_val, max_val], color='red', linestyle='--', label='Perfect Prediction')
    
        axes[i].set_xlabel('True Values')
        axes[i].set_ylabel('Predicted Values')
        axes[i+3].set_xlabel('True Values')
        axes[i+3].set_ylabel('Errors')
        axes[i].set_title(titles[i], y=titleposition)
        axes[i+3].set_title(titles[i+3], y=titleposition)
        axes[i].legend()
        axes[i].grid()
        axes[i+3].grid()
        axes[i].axis('equal')
        
    fig.savefig('./figures/visualize_results.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

if __name__ == "__main__":
    file_path1 = 'results/reprod.txt'
    file_path2 = 'results/errors.txt'
    inputs, true_values, predicted_values, errors = read_data(file_path1, file_path2)
    CRs = Rcoeff(true_values, predicted_values)
    print(f'Correlation Coefficients: {CRs}')
    visualize_results(inputs, true_values, predicted_values, errors, CRs)