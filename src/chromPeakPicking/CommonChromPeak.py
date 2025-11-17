import numpy as np


class ChromPeak:
    def __init__(
        self,
        peakApexIndex: int,
        peakApexRt: float,
        peakLeftFlankIndex: int,
        peakLeftFlankRt: float,
        peakRightFlankIndex: int,
        peakRightFlankRt: float,
        peakArea: float = None,
        peakApexHeight: float = None,
        peakFWHM: float = None,
        peakSNR: float = None,
        fileName: str = None,
        mz: float = None,
        eic: np.ndarray = None,
    ):
        self.peakApexIndex = peakApexIndex
        self.peakApexRt = peakApexRt
        self.peakLeftFlankIndex = peakLeftFlankIndex
        self.peakLeftFlankRt = peakLeftFlankRt
        self.peakRightFlankIndex = peakRightFlankIndex
        self.peakRightFlankRt = peakRightFlankRt
        self.peakArea = peakArea
        self.peakApexHeight = peakApexHeight
        self.peakFWHM = peakFWHM
        self.peakSNR = peakSNR
        self.fileName = fileName
        self.mz = mz
        if eic is not None:
            self.eic = eic  # Ensure eic is a numpy array of shape [2, :]
        else:
            self.eic = None

    @classmethod
    def from_indices_and_eic(cls, peakApexIndex: int, peakLeftFlankIndex: int, peakRightFlankIndex: int, eic: list, fileName: str = None, mz: float = None):
        # Extract retention times and intensities
        retention_times = eic[0, :]
        intensities = eic[1, :]

        # Calculate apex retention time and height
        peakApexRt = retention_times[peakApexIndex]
        peakApexHeight = intensities[peakApexIndex]

        # Calculate left and right flank retention times
        peakLeftFlankRt = retention_times[peakLeftFlankIndex]
        peakRightFlankRt = retention_times[peakRightFlankIndex]

        # Calculate baseline as the minimum intensity of the left and right flanks
        baseline = np.min(intensities[peakLeftFlankIndex:peakRightFlankIndex])

        # Calculate FWHM
        half_max = (peakApexHeight + baseline) / 2
        left_half_max_rt = None
        right_half_max_rt = None

        # Find the retention time at half max on the left side
        left_half_max_index = np.where(intensities[peakLeftFlankIndex : peakApexIndex + 1][::-1] <= half_max)[0]
        if left_half_max_index.size > 0:
            left_half_max_rt = retention_times[peakApexIndex - left_half_max_index[0]]
        else:
            left_half_max_rt = None

        # Find the retention time at half max on the right side
        right_half_max_index = np.where(intensities[peakApexIndex : peakRightFlankIndex + 1] <= half_max)[0]
        if right_half_max_index.size > 0:
            right_half_max_rt = retention_times[peakApexIndex + right_half_max_index[0]]
        else:
            right_half_max_rt = None

        # Calculate FWHM
        if left_half_max_rt is not None and right_half_max_rt is not None:
            peakFWHM = right_half_max_rt - left_half_max_rt
        else:
            peakFWHM = 0.0  # Default to 0 if FWHM cannot be calculated

        # Calculate peak area (simple trapezoidal integration)
        peakArea = sum((intensities[i] + intensities[i + 1]) / 2 * (retention_times[i + 1] - retention_times[i]) for i in range(peakLeftFlankIndex, peakRightFlankIndex))

        # Calculate SNR (signal-to-noise ratio)
        peakSNR = (peakApexHeight - baseline) / baseline if baseline > 0 else float("inf")

        return cls(peakApexIndex, peakApexRt, peakLeftFlankIndex, peakLeftFlankRt, peakRightFlankIndex, peakRightFlankRt, peakArea, peakApexHeight, peakFWHM, peakSNR, fileName, mz, eic)

    def __repr__(self):
        return (
            f"ChromPeak(peakApexIndex={self.peakApexIndex}, "
            f"peakApexRt={self.peakApexRt}, "
            f"peakLeftFlankIndex={self.peakLeftFlankIndex}, "
            f"peakLeftFlankRt={self.peakLeftFlankRt}, "
            f"peakRightFlankIndex={self.peakRightFlankIndex}, "
            f"peakRightFlankRt={self.peakRightFlankRt}, "
            f"peakArea={self.peakArea}, "
            f"peakApexHeight={self.peakApexHeight}, "
            f"peakFWHM={self.peakFWHM}, "
            f"peakSNR={self.peakSNR}, "
            f"fileName={self.fileName}, "
            f"mz={self.mz}, "
            f"eic={self.eic})"
        )


def getPeakStats(eic: np.ndarray, peakLeftFlankIndex: int, peakApexIndex: int, peakRightFlankIndex: int):
    peak = ChromPeak.from_indices_and_eic(peakApexIndex, peakLeftFlankIndex, peakRightFlankIndex, eic)
    return peak.peakFWHM, peak.peakArea, peak.peakSNR


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    # Generate random retention times
    retention_times = np.linspace(0, 10, 1000)

    # Generate a Gaussian-like peak
    peak_center = 5
    peak_width = 0.5
    peak_height = 100
    gaussian_peak = peak_height * np.exp(-((retention_times - peak_center) ** 2) / (2 * peak_width**2))

    # Add random noise
    noise_level = 5
    noise = np.random.normal(0, noise_level, retention_times.shape)
    intensities = gaussian_peak + noise

    # Ensure no negative intensities
    intensities = np.clip(intensities, 0, None)

    # Define the EIC as a 2D array
    eic = np.vstack((retention_times, intensities))

    # Define peak indices
    peakApexIndex = np.argmax(intensities)
    peakLeftFlankIndex = np.where(retention_times < peak_center - peak_width * 3)[0][-1]
    peakRightFlankIndex = np.where(retention_times > peak_center + peak_width * 3)[0][0]

    # Create ChromPeak object
    peak = ChromPeak.from_indices_and_eic(peakApexIndex, peakLeftFlankIndex, peakRightFlankIndex, eic)

    # Print the calculated FWHM
    fwhm, area, snr = getPeakStats(eic, peakLeftFlankIndex, peakApexIndex, peakRightFlankIndex)
    print("Peak Statistics:")
    print(f"   - FWHM: {fwhm}")
    print(f"   - Area: {area}")
    print(f"   - SNR: {snr}")

    # Plot the peak
    plt.plot(retention_times, intensities, label="EIC with noise")
    plt.axvline(retention_times[peakApexIndex], color="red", linestyle="--", label="Peak Apex")
    plt.axvline(retention_times[peakLeftFlankIndex], color="green", linestyle="--", label="Left Flank")
    plt.axvline(retention_times[peakRightFlankIndex], color="blue", linestyle="--", label="Right Flank")
    if peak.peakFWHM > 0 and peak.eic is not None:
        left_half_max_rt = peak.peakApexRt - peak.peakFWHM / 2
        right_half_max_rt = peak.peakApexRt + peak.peakFWHM / 2
        plt.axvspan(left_half_max_rt, right_half_max_rt, color="yellow", alpha=0.3, label=f"FWHM: {peak.peakFWHM:.2f}")
    plt.title("Extracted Ion Chromatogram (EIC)")
    plt.xlabel("Retention Time")
    plt.ylabel("Intensity")
    plt.legend()
    plt.show()
