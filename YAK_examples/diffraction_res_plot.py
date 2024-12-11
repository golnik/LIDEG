import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
from matplotlib.patches import FancyBboxPatch
import numpy as np
import scipy.integrate as integrate

# Constants
au2nm = 0.052917829614246
au2A = au2nm * 10
au2eV = 27.211396641308
au2Vnm = 5.14220826 * 10**2
au2fs = 0.02418884254

plt.rcParams.update({"font.size": 14})
plt.rcParams["mathtext.fontset"] = "cm"
fontdict = {"fontsize": 10}

# Load data
tfile = pd.read_csv("output/tfile.dat", sep="\s+", header=None, comment="#")
QM_TB_diffr_intra = pd.read_csv(
    "output/diffr_total_QM.dat", sep="\s+", header=None, comment="#"
)
TB_diffr_intra = pd.read_csv(
    "output/diffr_total_SC.dat", sep="\s+", header=None, comment="#"
)
YA_diffr_intra = pd.read_csv(
    "output/Yakovlev_data.dat", sep="\s+", header=None, comment="#"
)

x = tfile.iloc[:, 0]
x1 = YA_diffr_intra.iloc[:, 0] - YA_diffr_intra.iloc[:, 0][0] + 0.8
x2 = YA_diffr_intra.iloc[:, 0] - YA_diffr_intra.iloc[:, 0][0] + 0.2

# Data columns
Afield = tfile.iloc[:, 2] * au2Vnm / au2eV
Efield = tfile.iloc[:, 4] * au2Vnm
TB_D10 = TB_diffr_intra.iloc[:, 6]
TB_D11 = TB_diffr_intra.iloc[:, 7]
QM_TB_D10 = QM_TB_diffr_intra.iloc[:, 1]
QM_TB_D11 = QM_TB_diffr_intra.iloc[:, 2]
YA_D10 = YA_diffr_intra.iloc[:, 1]
YA_D11 = YA_diffr_intra.iloc[:, 4]


#########################

# Constants
sigma = 0.7213

# Make sure x and TB_D10 are numpy arrays and validate their shapes
assert len(x) == len(TB_D10), "x and TB_D10 must have the same length"

# Define the Gaussian kernel function
def gaussian_kernel(x, t, sigma):
    return np.exp(-((x - t) ** 2) / (sigma ** 2))

# Compute the transformed `New_TB_D10` using Gaussian convolution over limited range
New_TB_D10 = np.zeros_like(TB_D10)
New_TB_D11 = np.zeros_like(TB_D11)

# Number of points to look ahead and behind
num_points = 15

# Loop over each value of x and calculate the weighted sum within the range of previous and next 15 points
for i in range(len(x)):
    # Define the range of indices to consider
    start_idx = max(0, i - num_points)
    end_idx = min(len(x), i + num_points + 1)

    # Get the subarray of x and TB_D10 within the defined range
    sub_x = x[start_idx:end_idx]
    sub_TB_D10 = TB_D10[start_idx:end_idx]
    sub_TB_D11 = TB_D11[start_idx:end_idx]

    # If there are fewer points than expected, pad the subarray with the last available value
    if len(sub_x) < 2 * num_points + 1:
        pad_length = 2 * num_points + 1 - len(sub_x)
        sub_x = np.pad(sub_x, (0, pad_length), mode='edge')
        sub_TB_D10 = np.pad(sub_TB_D10, (0, pad_length), mode='edge')
        sub_TB_D11 = np.pad(sub_TB_D11, (0, pad_length), mode='edge')

    # Calculate Gaussian weights for each point in the subarray relative to x[i]
    weights = gaussian_kernel(x[i], sub_x, sigma)

    # Perform the weighted sum
    New_TB_D10[i] = np.sum(weights * sub_TB_D10)
    New_TB_D11[i] = np.sum(weights * sub_TB_D11)

# Replace the first 10 values of `TB_D10` with the computed `New_TB_D10` values
New_TB_D10[:10] = New_TB_D10[10]
New_TB_D11[:10] = New_TB_D11[10]


# Normalize data
TB_D10_filtted = New_TB_D10 / New_TB_D10[0] - 1
TB_D11_filtted = New_TB_D11 / New_TB_D11[0] - 1
# TB_D10_filtted = TB_D10 / TB_D10[0] - 1
# TB_D11_filtted = TB_D11 / TB_D11[0] - 1
QM_TB_D10_filtted = QM_TB_D10 / QM_TB_D10[0] - 1
QM_TB_D11_filtted = QM_TB_D11 / QM_TB_D11[0] - 1
YA_D10_filtted = YA_D10 / YA_D10[0] - 1
YA_D11_filtted = YA_D11 / YA_D11[0] - 1



#########################

x = x -15
x1 = x1 - 15
x2 = x2 - 15

# Create figure with extra space on the right for legends
fig, ax = plt.subplots(3, 1, figsize=(9, 6), sharex=True)

# Main plots
# Vector potential
ax[0].plot(x, Efield, color="tab:red", linestyle="--", label=r"Electric" "\n" r"Field")
ax[0].plot(x, Afield, color="tab:blue", label=r"Vector" "\n" r"Potential")
ax0 = ax[0].twinx()
ax0.plot(x, Afield,color="tab:blue", label=r"Vector" "\n" r"Potential")
ax[0].set_ylabel(r"Electric Field" "\n" r"[V/nm]")
ax0.set_ylabel(r"Vector Potential" "\n" r"[V$\cdot$ fs/m]")
ax[0].yaxis.set_major_locator(MaxNLocator(3))

ax0.set_yticks([])
ax0.set_ylim(-8, 8)
ax[0].set_ylim(-8, 8)

# Scattering intensity at [10]
ax[1].plot(x, TB_D10_filtted * 1e2, color="tab:red", label="SC-TRDI")
ax[1].plot(x, QM_TB_D10_filtted * 1e2, color="tab:purple",linestyle="-.", label="QM-TRDI")
ax[1].plot(x1, YA_D10_filtted * 1e2, color="tab:green", linestyle="--", label=r"Yakovlev" "\n" r"$\it{et}\ \it{al.}$ [15]")
ax[1].set_ylabel(r"Scattering" "\n" r"Intensity at [1,0]")
ax[1].text(
    -0.070,
    0.95,
    r"1E-2",
    transform=ax[1].transAxes,
    ha="left",
    va="bottom",
    fontdict=fontdict,
)
ax[1].yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
ax[1].yaxis.set_major_locator(MaxNLocator(3))
ax[1].set_ylim(-4.5, 4.5)

# Scattering intensity at [11]
ax[2].plot(x, TB_D11_filtted * 1e2, color="tab:red", label="SC-TRDI")
ax[2].plot(x, QM_TB_D11_filtted * 1e2, color="tab:purple",linestyle="-.", label="QM-TRDI")
ax[2].plot(x2, YA_D11_filtted * 1e2, color="tab:green", linestyle="--", label=r"Yakovlev" "\n" r"$\it{et}\ \it{al.}$ [15]")
ax[2].set_ylabel(r"Scattering" "\n" r"Intensity at [1,1]")
ax[2].text(
    -0.070,
    0.95,
    r"1E-2",
    transform=ax[2].transAxes,
    ha="left",
    va="bottom",
    fontdict=fontdict,
)
ax[2].yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
ax[2].yaxis.set_major_locator(MaxNLocator(3))
ax[2].set_ylim(-2.5, 2.5)
ax[2].set_xlim(0-15, 15)

# Add legends with a bordered box for each relevant subplot
for i, axi in enumerate(ax):
    if i == 0:
        # Add the legend for the first subplot only
        legend = axi.legend(
            loc="center left",
            bbox_to_anchor=(1.07, 0.5),
            frameon=False,
            handlelength=2,
            labelspacing=0.5,
        )

        # Create a surrounding FancyBboxPatch for the legend
        bbox = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        fancy_box = FancyBboxPatch(
            (bbox.x0, bbox.y0),
            bbox.width,
            bbox.height,
            boxstyle="round,pad=0.3",
            edgecolor="black",  # Black boundary
            facecolor="white",
            alpha=0.9,
            linewidth=1.5,  # Thickness of the boundary
            transform=fig.transFigure,
        )
        fig.patches.append(fancy_box)

    elif i == 2:
        # Add the legend for the third subplot and place it in the middle of the second and third plots
        legend = axi.legend(
            loc="center left",
            bbox_to_anchor=(1, 1.2),
            frameon=False,
            handlelength=2,
            labelspacing=0.5,
        )

        # # Create a surrounding FancyBboxPatch for the legend
        # bbox = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        # fancy_box = FancyBboxPatch(
        #     (bbox.x0, bbox.y0),
        #     bbox.width,
        #     bbox.height,
        #     boxstyle="round,pad=0.3",
        #     edgecolor="black",  # Black boundary
        #     facecolor="white",
        #     alpha=0.9,
        #     linewidth=1.5,  # Thickness of the boundary
        #     transform=fig.transFigure,
        # )
        # fig.patches.append(fancy_box)

# Final adjustments
fig.subplots_adjust(left=0.11, right=0.775, bottom=0.10, top=0.95, hspace=0.3)
fig.text(0.5, 0.01, "Time [fs]", ha="center")
plt.savefig("figures/diffr_compare_total_with_boxed_legends.pdf", format="pdf")
