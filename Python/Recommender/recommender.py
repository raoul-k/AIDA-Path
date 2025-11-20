
import pandas as pd

# Arguments:
# df_list ... list of pandas dataframes (binary matrices of disorders)
# has_symptoms ... list of symptoms(strings)
def recommend_next_symptoms(df_list, has_symptoms):
    df_all_disorders = pd.concat(df_list)
    df_all_disorders = df_all_disorders.fillna(0)

    df_filtered = df_all_disorders
    string_query = ' and '.join([s+'==1' for s in has_symptoms])
    df_filtered = df_filtered.query(string_query)
    df_filtered = df_filtered.drop(has_symptoms, axis=1)
    disorders_coulds = set(df_filtered['disorder'].tolist())
    

    #Output 1: Possible Disorders
    print("Possible disorders:",", ".join(disorders_coulds))
    print()
    
    df_grouped  = df_filtered.groupby('disorder', as_index=False).agg('mean')
    df_grouped['1'] = df_grouped.apply(lambda row: row[row == 1].index.tolist(), axis=1)
    df_grouped['0'] = df_grouped.apply(lambda row: row[row == 0].index.tolist(), axis=1)
    
    #alternative: set(itertools.chain.from_iterable(df_grouped['1'].tolist()))
    diff1 = set(sum(df_grouped['1'].tolist(),[]))
    diff0 = set(sum(df_grouped['0'].tolist(),[]))
    symptoms2check = diff0 & diff1
    
    #Output 2: Important Symptoms, (One-Sided), Two-Sided
    #print("Important symptoms to check for presence:",", ".join(diff1))
    #print("Important symptoms to check for absence:",", ".join(diff0))
    print("Important symptoms to check:",", ".join(symptoms2check))
    print()
    
    for s in symptoms2check:
        disorders0 = set()
        disorders1 = set()
        for index, row in df_grouped.iterrows():
            if s in row['1']:
                disorders1.add(row['disorder'])
            if s in row['0']:
                disorders0.add(row['disorder'])
        #Output 3: Symptom/Disorder Pairings
        print("Check", s, "to differentiate between", ", ".join(disorders1), "and", ", ".join(disorders0))