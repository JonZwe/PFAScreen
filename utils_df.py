
import pandas as pd
import pyopenms as oms
import matplotlib.pyplot as plt
import pyperclip

"""
Collection of various functions applicable to the pandas DataFrame objects used in PFAScreen.
Jonathan Zweigle, 06/2025
"""

def get_by_accurate_mass(df, mz, tol):
    """
    Function to get all entries in a dataframe within a certain m/z tolerance
    """
    if isinstance(mz, str):
        mz = oms.EmpiricalFormula(mz).getMonoWeight()

    condition = (df['mz'] >= mz - tol) & (df['mz'] <= mz + tol)
    return df[condition]


def filter_df(df, column, lower=None, upper=None, reset_index=True, keep_nan=False):
    """Filter a pandas DataFrame by column range."""
    condition = pd.Series(True, index=df.index)

    if lower is not None:
        condition &= df[column] >= lower
    if upper is not None:
        condition &= df[column] <= upper

    if keep_nan:
        condition |= df[column].isna()

    print(f'{len(df)} -> {len(df[condition])}')

    if reset_index:
        return df[condition].reset_index(drop=True)
    return df[condition]


def get_n_highest(df, col, n=20):
    return df.sort_values(by=col, ascending=False).iloc[:n]


def copy_to_clipboard(df, col):
    pyperclip.copy(", ".join(df[col].astype(str).to_list()))


def fold_change_filter_simple(df, sample_names, blank, fold_change):

    """
    Drop rows where the fold change of any sample compared to the blank is below a threshold.
    """

    l = len(df)
    idx = df.index[(df[sample_names].drop(columns=[blank]).div(df[sample_names][blank], axis=0) >= fold_change).any(axis=1)]
    
    df = df.loc[idx].reset_index(drop=True)
    print(f'{l} -> {len(df)}')
    return df


def fold_change_filter(df, sample_names, blank, fold_change,
                       blank_aggregation='mean'):
    """
    Keep rows where at least one non-blank sample is at least `fold_change`
    times the aggregated abundance across the blank samples.

    Parameters
    ----------
    blank : str or list[str]
        One blank sample name or a list of blank sample names.
    blank_aggregation : {'mean', 'median', 'max'}
        Method used to combine multiple blank abundances.
    """
    blank_samples = [blank] if isinstance(blank, str) else list(blank)
    non_blank_samples = [
        sample for sample in sample_names if sample not in blank_samples
    ]

    if not non_blank_samples:
        raise ValueError("sample_names must contain at least one non-blank sample.")

    missing_samples = set(blank_samples + non_blank_samples) - set(df.columns)
    if missing_samples:
        raise KeyError(f"Samples not found in DataFrame: {sorted(missing_samples)}")

    if blank_aggregation == 'mean':
        blank_reference = df[blank_samples].mean(axis=1)
    elif blank_aggregation == 'median':
        blank_reference = df[blank_samples].median(axis=1)
    elif blank_aggregation == 'max':
        blank_reference = df[blank_samples].max(axis=1)
    else:
        raise ValueError(
            "blank_aggregation must be 'mean', 'median', or 'max'."
        )

    fold_changes = df[non_blank_samples].div(blank_reference, axis=0)
    keep_mask = fold_changes.ge(fold_change).any(axis=1)

    initial_count = len(df)
    df = df.loc[keep_mask].reset_index(drop=True)

    print(f'{initial_count} -> {len(df)}')
    return df


def plot_feature_corr(
    df,
    feature_1,
    feature_2,
    sample_names,
    identifier='id',
    index_1=0,
    index_2=0,
    tol=0.005,
    annotate=False
):
    """
    Plot sample intensities for two features.

    identifier:
        'id' or 'mz'
    """

    sample_names = list(dict.fromkeys(sample_names))

    if identifier == 'id':
        fid_1 = feature_1
        fid_2 = feature_2

    elif identifier == 'mz':
        hits_1 = get_by_accurate_mass(df, feature_1, tol)
        hits_2 = get_by_accurate_mass(df, feature_2, tol)

        if hits_1.empty:
            raise ValueError(f'No feature found near m/z {feature_1}.')
        if hits_2.empty:
            raise ValueError(f'No feature found near m/z {feature_2}.')

        if len(hits_1) > 1:
            print(f'More than one feature found for m/z: {feature_1}')
        if len(hits_2) > 1:
            print(f'More than one feature found for m/z: {feature_2}')

        fid_1 = hits_1.index[index_1]
        fid_2 = hits_2.index[index_2]

    else:
        raise ValueError("identifier must be either 'id' or 'mz'.")

    x = df.loc[fid_1, sample_names]
    y = df.loc[fid_2, sample_names]

    # Keep only samples where both features were detected.
    valid = x.notna() & y.notna()
    x_valid = x[valid]
    y_valid = y[valid]

    if len(x_valid) < 2:
        raise ValueError(
            'At least two shared non-missing samples are required for plotting.'
        )

    correlation = x_valid.corr(y_valid)

    plt.figure(figsize=(7, 6))
    plt.scatter(x_valid, y_valid)

    if annotate:
        for sample, x_value, y_value in zip(
            x_valid.index,
            x_valid.values,
            y_valid.values
        ):
            plt.annotate(
                sample,
                (x_value, y_value),
                xytext=(4, 4),
                textcoords='offset points'
            )

    mz_1 = df.loc[fid_1, 'mz']
    mz_2 = df.loc[fid_2, 'mz']

    plt.xlabel(f'Feature {fid_1} | m/z {mz_1:.5f}')
    plt.ylabel(f'Feature {fid_2} | m/z {mz_2:.5f}')
    plt.title(
        f'Correlation: r = {correlation:.3f}, '
        f'n = {len(x_valid)} shared samples'
    )
    plt.tight_layout()
    plt.show()

    return {
        'feature_id_1': fid_1,
        'feature_id_2': fid_2,
        'mz_1': mz_1,
        'mz_2': mz_2,
        'correlation': correlation,
        'n_shared': len(x_valid),
        'shared_samples': list(x_valid.index)
    }


def find_correlating_features(df, sample_names, mz, threshold=0.8, index=0, tol=0.005):

    sample_names = list(dict.fromkeys(sample_names))

    i = get_by_accurate_mass(df, mz, tol=tol)
    if len(i) > 1:
        print('More than one feature found for mz:', mz)
    i = i.index[index]

    ints = df.loc[i, sample_names]
    ref_compound = df.loc[i, 'compound_names'] if 'compound_names' in df.columns else None

    results = []
    for j in range(len(df)):
        ints_j = df.iloc[j][sample_names]
        n = (ints.notna() & ints_j.notna()).sum()
        corr = ints.corr(ints_j)
        if corr > threshold and n > 2:
            results.append({
                'feature_idx': df.index[j],
                'mz': df.iloc[j]['mz'],
                'rt': df.iloc[j]['rt'],
                'corr': corr,
                'n': n
            })

    result_df = pd.DataFrame(results).sort_values('corr', ascending=False)
    return result_df


def filter_alignment_table(df_data, df_meta,
                           group_col,
                           sample_col,
                           blank_label,
                           fold_change_threshold=5.0,
                           fold_change_type='mean',
                           rel_replicate_abundance=0.1,
                           max_replicate_rsd=0.75,
                           set_zeros_to_nan=True):
    """

    # NOTE: This function still needs to be rigerously tested and validated (04/08/2024)!

    Filters an alignment table based on fold change (vs blanks), relative abundance, and replicate RSD.
    Failing values are masked (set to NaN) on a per-group basis.
    
    Parameters:
        df_data (DataFrame): intensity matrix (rows = features, cols = samples)
        df_meta (DataFrame): metadata with sample group assignments
        rel_replicate_abundance (float): min proportion of max intensity within group
        max_replicate_rsd (float): max allowed RSD across replicates
        fold_change_threshold (float): min fold change over blanks
        fold_change_type (str): 'mean' or 'max' (blank reference)
        group_col (str): column in df_meta indicating group
        sample_col (str): column in df_meta matching df_data columns
        blank_label (str): label in group_col for blank samples
        drop_all_nan_rows (bool): drop features that are fully masked (NaN)
    
    Returns:
        df_filtered (DataFrame): masked intensity matrix
        final_mask (DataFrame): boolean mask of kept values
    """

    print(f'Initial number of features: {len(df_data)}')

    # Initialize full mask (all True)
    full_mask = pd.DataFrame(True, index=df_data.index, columns=df_data.columns)
    total_features = len(full_mask)

    # 1. Fold change filtering — apply first!
    if fold_change_threshold:

        blank_cols = df_meta[df_meta[group_col] == blank_label][sample_col].values
        non_blank_groups = df_meta[df_meta[group_col] != blank_label][group_col].unique()

        fc_mask = pd.DataFrame(True, index=df_data.index, columns=df_data.columns)

        for group in non_blank_groups:
            group_cols = df_meta[df_meta[group_col] == group][sample_col].values
            sample_mean = df_data[group_cols].mean(axis=1)

            if fold_change_type == 'mean':
                blank_mean = df_data[blank_cols].mean(axis=1).replace(0, np.nan)
                fc = sample_mean / blank_mean
            elif fold_change_type == 'max':
                blank_max = df_data[blank_cols].max(axis=1).replace(0, np.nan)
                fc = sample_mean / blank_max
            else:
                raise ValueError("fold_change_type must be 'mean' or 'max'")

            # Mask all values in group if FC fails
            fail_mask = fc < fold_change_threshold
            fc_mask.loc[fail_mask, group_cols] = False

        full_mask &= fc_mask
        print('Fold change filtering \n')
        
        [print(f"{group}: {total_features - (~full_mask[df_meta[df_meta[group_col] == group][sample_col].values]).all(axis=1).sum()}") 
        for group in df_meta[group_col].unique()]
    else:
        print('Fold change filtering skipped.')

    # 2. Relative abundance filtering
    if rel_replicate_abundance:
        for group in df_meta[group_col].unique():
            group_cols = df_meta[df_meta[group_col] == group][sample_col].values
            sub = df_data[group_cols]
            max_vals = sub.max(axis=1)
            threshold_vals = (max_vals * rel_replicate_abundance).values[:, np.newaxis]

            group_mask = sub >= threshold_vals
            any_fail = ~group_mask.all(axis=1)
            group_mask.loc[any_fail] = False

            full_mask.loc[:, group_cols] &= group_mask

        print('Relative replicate abundance \n')
        [print(f"{group}: {total_features - (~full_mask[df_meta[df_meta[group_col] == group][sample_col].values]).all(axis=1).sum()}") 
        for group in df_meta[group_col].unique()]
    else:
        print('Relative replicate abundance filtering skipped.')

    # 3. RSD filtering
    if max_replicate_rsd:
        for group in df_meta[group_col].unique():
            group_cols = df_meta[df_meta[group_col] == group][sample_col].values
            sub = df_data[group_cols]
            mean = sub.mean(axis=1)
            std = sub.std(axis=1)
            rsd = (std / mean).replace([np.inf, -np.inf], np.nan)

            fail_mask = rsd > max_replicate_rsd
            full_mask.loc[fail_mask, group_cols] = False

        print('Replicate RSD filtering \n')
        [print(f"{group}: {total_features - (~full_mask[df_meta[df_meta[group_col] == group][sample_col].values]).all(axis=1).sum()}") 
        for group in df_meta[group_col].unique()]
    else:
        print('Replicate RSD filtering skipped.')

    # 4. Apply full mask
    df_filtered = df_data.where(full_mask)
    
    # 5. Set zeros to NaN
    if set_zeros_to_nan:
        df_filtered = df_filtered.replace([0, 0.0], np.nan)

    print(f'Final number of features retained: {len(df_filtered)}')

    return df_filtered, full_mask