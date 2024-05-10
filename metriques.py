import csv

all_genus = ['Achromobacter', 'Acinetobacter', 'Actinobacillus', 'Aerococcus', 'Aeromonas', 'Aliarcobacter', 'Burkholderia', 'Campylobacter', 'Citrobacter', 'Clostridioides', 'Corynebacterium', 'Desulfovibrio', 'Enterobacter', 'Enterococcus', 'Escherichia', 'Haemophilus', 'Helicobacter', 'Herbaspirillum', 'Klebsiella', 'Kluyvera', 'Lachnoclostridium', 'Leclercia', 'Lelliottia', 'Listeria', 'Morganella', 'Mycobacterium', 'Mycolicibacillus', 'Mycolicibacterium', 'Neisseria', 'Phytobacter', 'Pluralibacter', 'Proteus', 'Providencia', 'Pseudomonas', 'Raoultella', 'Salmonella', 'Serratia', 'Shigella', 'Staphylococcus', 'Streptococcus', 'Stutzerimonas', 'Vibrio', 'Yersinia']
all_db = ["argannot", "card", "ncbi", "resfinder"]

l_data_2 = []
for i in all_db:
    l_data = []
    data_f = {
        i : ""
    }
    for y in all_genus:
        data = {
            "NC": 0,
            "TP": 0,
            "TN": 0,
            "FN": 0,
            "FP": 0,
        }
        with open(f"./data/genome/{i}/{y}.txt", "r") as file:
            reader = csv.DictReader(file, delimiter='\t')
            raw_data = [line for line in reader]

            for z in raw_data:
                 if z["relax"] == "NC":
                     data["NC"] += 1
                 if z["relax"] == "TP":
                     data["TP"] += 1
                 if z["relax"] == "TN":
                     data["TN"] += 1
                 if z["relax"] == "FP":
                     data["FP"] += 1
                 if z["relax"] == "FN":
                     data["FN"] += 1

        l_data.append(data)
    data_f[i] = l_data
    l_data_2.append(data_f)

l_metrique_f = []
for i in l_data_2:
    l_metrique = []
    for key, value in i.items():
        data = value  # Récupérer les données correspondant à la clé actuelle
        for y in data:
            metrique = {
                "sensibilité": [],
                "specificity": [],
                "accuracy": [],
                "negative predict": [],
                "precision": []
            }
            fp = y["FP"]
            fn = y["FN"]
            tp = y["TP"]
            tn = y["TN"]

            if tp + fn != 0:
                sensi = tp / (tp + fn)
            else:
                sensi = 0
            if tn + fp != 0:
                speci = tn / (tn + fp)
            else:
                speci = 0
            if tp + tn + fn + fp != 0:
                accuracy = (tp + tn) / (tp + tn + fn + fp)
            else:
                accuracy = 0
            if tn + fn != 0:
                neg = tn / (tn + fn)
            else:
                neg = 0
            if tp + fp != 0:
                precision = tp / (tp + fp)
            else:
                precision = 0

            metrique["specificity"].append(speci)
            metrique["sensibilité"].append(sensi)
            metrique["accuracy"].append(accuracy)
            metrique["negative predict"].append(neg)
            metrique["precision"].append(precision)
            l_metrique.append(metrique)
        l_metrique_f.append(l_metrique)

for i in l_metrique_f:
    l_sensi = []
    l_speci = []
    l_accu = []
    l_neg = []
    l_precision = []
    for j in i:
        l_sensi.append(j["sensibilité"])
        l_speci.append(j["specificity"])
        l_accu.append(j["accuracy"])
        l_neg.append(j["negative predict"])
        l_precision.append(j["precision"])
    l_sensi = [item for sublist in l_sensi for item in sublist]
    l_speci = [item for sublist in l_speci for item in sublist]
    l_accu = [item for sublist in l_accu for item in sublist]
    l_neg = [item for sublist in l_neg for item in sublist]
    l_precision = [item for sublist in l_precision for item in sublist]
    print(f"sensi {l_sensi}")
    print(f"speci {l_speci}")
    print(f"accu {l_accu}")
    print(f"neg {l_neg}")
    print(f"preci {l_precision}")
    print("==============")

