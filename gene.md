## Gene map

Read the NBCI downloaded transcript with biopython

    from Bio import SeqIO

    transcript = SeqIO.read('NM_004269.4.gb', 'genbank')

Convert it to a pandas dataframe, namely
in a `Bio.SeqRecord.SeqRecord`, features are stored as a list of `Bio.SeqFeature.SeqFeature` in `record.features`.
The `qualifiers` dictionary within the feature has always a list, so those will be lists.

    import Bio
    import pandas as pd
    
    def convert(feat:Bio.SeqFeature.SeqFeature) -> dict:
        return {'id': feat.id,
                'type': feat.type,
                'location_strand': feat.location.strand, 
                'location_start': feat.location.start,
                'location_end': feat.location.end,
                **feat.qualifiers}
                
    
    feats = pd.DataFrame([convert(feat) for feat in transcript.features])

Plotting with plotly, whereas exons are to scale, spacer is a fake length for a half intron.

    import plotly.graph_objects as go
    
    spacer = 100
    
    exo_feats = feats.loc[feats.type == 'exon'].reset_index()
    spacing = exo_feats.index.to_series() * spacer * 2 + spacer
    
    get_zeros = lambda table: pd.Series([0] * len(table))
    
    fig = go.Figure(data=[go.Candlestick(x=get_zeros(exo_feats),
                                        open=exo_feats.location_start + spacing,
                                        high=exo_feats.location_end + spacing + spacer,
                                        low=exo_feats.location_start + spacing - spacer,
                                        close=exo_feats.location_end + spacing, 
                    # orientation='h' this does not work.
                                        ),
                         ])
    fig.update_layout(xaxis_rangeslider_visible=False)
    fig.show()
    
This was modified to included

    go.Scatter(x=get_zeros(variants),
               y=variants.position + 23,
               mode='markers')
               
It was saved `fig.write_image("candlesticks2.svg")` and rotated in Illustrator.