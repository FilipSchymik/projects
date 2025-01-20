def Fischer_enrich(background_df, sample_df):
    '''
    Fischer's exact test for go-term enrichment analysis.
    Input:
        background_terms: dataframe with go-term IDs, names, and counts in the background set
        target_terms: dataframe with go-term IDs, names, and counts in the target set
    '''

    results = []

    # Get unique proteins in each set
    background_proteins = set(background_df["protein"].unique())
    sample_proteins = set(sample_df["protein"].unique())

    # Get unique GO terms for each aspect
    all_go_terms = set(background_df["go_term_id"].unique())

    for aspect in ["BP", "MF", "CC"]:
        for go_term in all_go_terms:
            # Check if the GO term is present in the current aspect
            if not ((background_df['go_term_id'] == go_term) & (background_df['aspect'] == aspect)).any():
                continue