# Removing 50 Hz mains interference from an ECG with an FIR notch filter

BEC504 Digital Signal Processing - Experiential Learning, Question 13
V Semester, A Section, ECE, NIE Mysuru, AY 2026-27
Course instructor: Dr. Ashok K

Team 2
- Gagan Manjunath (4NI24EC047)
- Manoj Girish Bhat (4NI24EC071)
- B Srivatsa Joshi (4NI24EC026)

## About
A raw ECG picks up 50 Hz hum from the mains. This project removes it with a linear-phase FIR notch filter in MATLAB and checks how much of the ECG is disturbed in the process.

The filter is a 201-tap (N = 200) Hamming-window band-stop filter with a stop band of 45-55 Hz, designed from the window-method formula and checked against `fir1`. The FIR difference equation is also written out with loops and checked against `filter()`.

The program covers the five parts of the question:
1. design of the notch filter
2. magnitude response near 50 Hz
3. 50 Hz component before and after filtering
4. attenuation at 50 Hz and distortion of the ECG
5. effect of filter order (50-400) and notch width (2-16 Hz)

## Data
ECG-ID Database from PhysioNet, Person_25, record rec_1 (raw ECG, lead I), sampled at 500 Hz for 20 s.
https://physionet.org/content/ecgiddb/1.0.0/Person_25/

The raw channel has a clear mains tone at about 50.35 Hz. The data files are in `data/`.

## How to run
Open MATLAB in this folder (Signal Processing Toolbox needed) and run

    q13_fir_notch

The figures, tables and the console output are saved in `results/`.

## Results
| | |
|---|---|
| Attenuation at 50 Hz | 48.1 dB |
| Reduction of the 50 Hz peak in the recording | 49.9 dB |
| Gain over 0.5-40 Hz | within 0.03 dB |
| R-peak amplitude change | +0.49 % |
| Normalised RMS error against the ideal-cleaned signal | 2.91 % |
| Delay | 100 samples (200 ms) |

The 100 Hz harmonic of the mains is not removed by this filter.

## Files
- `q13_fir_notch.m` - MATLAB program
- `data/` - ECG recording (.dat and .hea)
- `results/` - plots (fig1 to fig6), comparison tables (.csv) and console output
- `Q13_FIR_Notch_50Hz_Report.pdf` and `.docx` - the report
