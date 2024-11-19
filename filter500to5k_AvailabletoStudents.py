from scipy.signal import butter

def filter500to5k(Fs):
    """
    Returns a discrete time filter coefficients in sos format
    :param Fs: sampling frequency of data
    :return: sos, G
    """

    N = 1 # order, for bandpass filter the effective order is doubled
    fc1 = 500 # cutoff frequency 1
    fc2 = 5000 # cutoff frequency 2

    # Design the Butterworth filter
    sos = butter(N, [fc1, fc2], btype='band', fs=Fs, output='sos')
    return sos