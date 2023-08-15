import csv
from scipy.signal import find_peaks, peak_widths, peak_prominences

def get_peaks(results_dict, original_record_dict, fasta_file_name, db_file_name):

    peaks_list = []

    for key, arr in results_dict.items():
        peaks, _ = find_peaks(arr, width=16)
        prom = peak_prominences(arr, peaks)
        rh = 0.65
        widths, _, _, _ = peak_widths(arr, peaks, rh, prom)
        if peaks.size != 0:
            for i, peak in enumerate(peaks):
                start_idx = int(peak - (widths[i] / 2) + 1)
                if start_idx < 0:
                    start_idx = 0
                end_idx = int(peak + (widths[i] / 2) + 1)  # Add 1 to include the last index
                peak_string = original_record_dict[key][start_idx:end_idx]

                total_length = end_idx - start_idx
                if total_length > 27:
                    sub_length = total_length - 27
                    end_idx = end_idx - (sub_length)

                if len(peak_string) > 16:
                    peaks_list.append([key, arr[peak], start_idx, end_idx, peak_string[:27]])
                else:
                    peaks_list.append([key, '0', '0', 'Zero Reads',])    
        else:
            peaks_list.append([key, '0', '0', 'Zero Reads'])

        peaks = []
        widths = []

    out_file = f"{fasta_file_name}({db_file_name}).csv"

    with open(out_file, 'w') as f:
        fields = ['miRNA ID', 'Expression Level','Peak Start', 'Peak End', 'Peak Data']
        writer = csv.writer(f)
        writer.writerow(fields)
        writer.writerow([])

        for i in peaks_list:
            writer.writerow(i)


