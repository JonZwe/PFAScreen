# %%
import numpy as np
import matplotlib.pyplot as plt
from pylab import figure

def MS2_extractor(
        Df_MS2RawData,
        Df_FindPFAS,
        prec_mz,
        ):
    
    if np.sum(np.abs(Df_MS2RawData['m/z'] - prec_mz) < 0.005) == 0: # NOTE: Make mass tol as user input
        print('No MS2 available')
    else:
        idx = np.abs(Df_MS2RawData['m/z'] - prec_mz) < 0.005

        prec_mz_true = Df_MS2RawData['m/z'][idx].values[0]
        mz_array = Df_MS2RawData['MS2SpecMZ'][idx].values[0]
        intens_array = Df_MS2RawData['MS2SpecIntens'][idx].values[0]

        intens_array_normalized = intens_array/np.max(intens_array)
        idx_label = intens_array_normalized > 0.05
        mz_array_label = mz_array[idx_label]
        intens_array_label = intens_array[idx_label]

        figure()
        _, sl1, _ = plt.stem(mz_array, intens_array, 'Black', markerfmt=" ", basefmt=" ")
        plt.setp(sl1, color = 'Black', linewidth = 1)
        for i, txt in enumerate(np.round(mz_array_label, 4)):
            plt.annotate(txt, (mz_array_label[i],intens_array_label[i]), color = 'Black', rotation = 0, fontsize = 9)

        # Find MS/MS spectrum in Df_FindPFAS
        if type(Df_FindPFAS) != str:
            idx_FindPFAS = Df_FindPFAS['m/z_MSMS'] == prec_mz_true
            if np.sum(idx_FindPFAS) > 0: # CHANGED FROM == 1! NEED TO BE VALIDATED!
                start_idx = Df_FindPFAS.columns.get_loc('formula_diagnostic') + 1 # find first index with diff (e.g. CF2)
                end_idx = Df_FindPFAS.columns.get_loc(f'mz_{Df_FindPFAS.columns[start_idx]}')

                diffs = list(Df_FindPFAS.columns[start_idx:end_idx].values) # get diffs
                for diff in diffs:
                    if len(Df_FindPFAS[idx_FindPFAS][f'mz_{diff}'].values[0]) > 0:
                        _, sl2, _ = plt.stem(Df_FindPFAS[idx_FindPFAS][f'mz_{diff}'].values[0], 
                                             Df_FindPFAS[idx_FindPFAS][f'intens_{diff}'].values[0], markerfmt=" ", basefmt=" ")
                        plt.setp(sl2, linewidth = 1.5)

                if len(Df_FindPFAS[idx_FindPFAS]['mz_peaks_diagnostic'].values[0]) > 0:
                    _, sl3, _ = plt.stem(Df_FindPFAS[idx_FindPFAS]['mz_peaks_diagnostic'].values[0], 
                                         Df_FindPFAS[idx_FindPFAS]['intens_peaks_diagnostic'].values[0], 'Red', markerfmt=" ", basefmt=" ")
                    plt.setp(sl3, color = 'Red', linewidth = 1.5)

                    label = [n +'\n' for n in Df_FindPFAS[idx_FindPFAS]['formula_diagnostic'].values[0].to_list()]
                    for i, txt in enumerate(label):
                        plt.annotate(txt, (Df_FindPFAS[idx_FindPFAS]['mz_peaks_diagnostic'].values[0][i],
                                           Df_FindPFAS[idx_FindPFAS]['intens_peaks_diagnostic'].values[0][i]),
                                           color = 'Red', rotation = 0, fontsize = 9)

                if isinstance(Df_FindPFAS['new_formulas'][idx_FindPFAS].values[0], list): # NOTE: There could be a bug here, check again
                    #label2 = [n +'\n\n' for n in Df_FindPFAS['new_formulas'][idx_FindPFAS].values[0]] # NOTE: include two linebreaks
                    u, c = np.unique(Df_FindPFAS['frag_idx'][idx_FindPFAS].values[0], return_counts = True)
                    prefix_c = np.zeros(len(Df_FindPFAS['frag_idx'][idx_FindPFAS].values[0]))
                    for x in range(len(u)):
                        prefix_c[np.where(u[x] == Df_FindPFAS['frag_idx'][idx_FindPFAS].values[0])[0]] = np.arange(0, c[x])
                    prefix_c = prefix_c + 2
                    prefix_c = prefix_c.astype('int')
                    label2 = []
                    for n, formula in enumerate(Df_FindPFAS['new_formulas'][idx_FindPFAS].values[0]):
                        label2.append(formula + prefix_c[n]*'\n')

                    for i, txt in enumerate(label2):
                        plt.annotate(txt, (Df_FindPFAS['mz_peaks'][idx_FindPFAS].values[0][Df_FindPFAS['frag_idx'][idx_FindPFAS].values[0]][i],
                                           Df_FindPFAS['intens_peaks'][idx_FindPFAS].values[0][Df_FindPFAS['frag_idx'][idx_FindPFAS].values[0]][i]),
                                           color = 'Blue', rotation = 0, fontsize = 9)

                    #str_formulas = ' '.join(Df_FindPFAS['new_formulas'][idx_FindPFAS].values[0])
                    #plt.title(f'MS2 spectrum (m/z = {np.round(prec_mz_true, 4)}) \n {str_formulas}')
                    plt.title(f'MS2 spectrum (m/z = {np.round(prec_mz_true, 4)})')

        else:
            plt.title(f'MS2 spectrum (m/z = {np.round(prec_mz_true, 4)})')
        
        plt.ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0), useMathText=True)
        plt.xlabel('m/z')
        plt.ylabel('Counts (-)')
        plt.ylim(ymin = 0)
        plt.show()


        mz = mz_array[mz_array < prec_mz_true + 5]
        intens = intens_array[mz_array < prec_mz_true + 5]
        MD = mz - np.round_(mz, decimals = 0)  # NOTE: Change rounding!!
        MD_prec = prec_mz_true - np.round_(prec_mz_true, decimals = 0)
        idx = intens > 100

        figure()
        cmap = 'Reds'
        plt.scatter(mz[idx], MD[idx], 80, c = np.log10(intens[idx]), cmap = cmap, alpha = 0.8, edgecolors='lightgrey')
        plt.axhline(y = MD_prec, color = 'grey', linestyle = ':', alpha = 0.5)
        cbar = plt.colorbar()
        cbar.set_label('log(Intensity)')
        plt.xlabel('m/z')
        plt.ylabel('MD')
        plt.title(np.round(prec_mz_true, 4))
        plt.show()