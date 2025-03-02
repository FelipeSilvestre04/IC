import sys
import os

sys.path.append(os.path.abspath("/home/fsilvestre"))
from Cutting_Stock_Problem.Ambiente.Main.CSP_embraer import CSP, ler_poligonos, suavizar_poligono, tratar_lista, ajustar_poligono
from Cutting_Stock_Problem.Algoritimos.Heuristicas.GCG.Main.GCG_embraer import GCG
from Cutting_Stock_Problem.Algoritimos.Metaheuristicas.GRASP.Main.GRASP_GCG_embraer import GRASP_GCG, area_fecho_retangular, porcentagem_area
from Cutting_Stock_Problem.Algoritimos.Metaheuristicas.GRASP.Main.GRASP_GCG_embraer import PackS_GG
from Cutting_Stock_Problem.Algoritimos.Metaheuristicas.GRASP.Main.GRASP_embraer import GRASP, calcular_fecho_area, calcular_area_preenchida
from Cutting_Stock_Problem.Algoritimos.Metaheuristicas.GRASP.Main.GRASP_embraer import PackS_G
from Cutting_Stock_Problem.Algoritimos.Metaheuristicas.RKGA.Main.RKGA import RKGA, pre_processar_NFP

from Cutting_Stock_Problem.Algoritimos.Metaheuristicas.BRKGA.Main.BRKGA import BRKGA
from Cutting_Stock_Problem.Algoritimos.Heuristicas.GreedyIniSol.Main.GreedyInisol import GreedyInisol
from Cutting_Stock_Problem.Algoritimos.Heuristicas.GreedySearch.Main.GreedySearch import GreedySearch

from datetime import datetime, date
import time
import copy
import csv

# Função para criar a pasta com base na data e hora atuais
def criar_pasta_teste():
    now = datetime.now()
    folder_name = f"testes_{now.strftime('%Y%m%d_%H%M%S')}"
    os.makedirs(folder_name, exist_ok=True)
    return folder_name

def escrever_poligonos(poligonos, arquivo_saida):
    with open(arquivo_saida, 'w') as f:
        f.write(f"{len(poligonos)}\n\n")
        for i, poligono in enumerate(poligonos):
            f.write(f"{len(poligono)}\n")
            for x, y in poligono:
                f.write(f"{int(x)} {int(y)}\n")
            if i != len(poligonos) - 1:
                f.write("\n")

def create_results_file(folder):
    current_time = datetime.now()
    filename = os.path.join(folder, f"algoritimos_finais_{current_time.day}_{current_time.month}_{current_time.hour}_{current_time.minute}.txt")
    with open(filename, "w", encoding='utf-8') as file:
        file.write("="*190 + "\n")
        file.write(f"RESULTADOS GERAIS: {current_time.strftime('%d/%m/%Y %H:%M')}\n")
        file.write("="*190 + "\n\n")
    return filename

def calcular_area(poligono):
    if poligono[0] != poligono[-1]:
        poligono = poligono + [poligono[0]]
    n = len(poligono)
    area = 0.0
    for i in range(n-1):
        x1, y1 = poligono[i]
        x2, y2 = poligono[i+1]
        area += (x1 * y2) - (x2 * y1)
    return abs(area) / 2.0

def configurar_teste():
    # Solicita quantas dimensões serão testadas (qualquer número maior que 0)
    while True:
        try:
            num_dimensoes = int(input("Digite o número de dimensões a serem testadas: "))
            if num_dimensoes > 0:
                break
            print("O número de dimensões deve ser maior que 0.")
        except ValueError:
            print("Por favor, digite um número válido.")
    
    # Para cada dimensão, o usuário informa a porcentagem (ex: 50 para 50%)
    dimensoes_selecionadas = []
    for i in range(num_dimensoes):
        while True:
            dim = input(f"Digite a {i+1}ª dimensão (em porcentagem, sem o símbolo %): ")
            try:
                valor = float(dim)
                if valor > 0:
                    dimensoes_selecionadas.append(valor)
                    break
                else:
                    print("Por favor, digite um valor maior que zero.")
            except ValueError:
                print("Por favor, digite um número válido.")
    
    # Seleção dos algoritmos
    print("\nSeleção de algoritmos (1 para executar, 0 para não executar):")
    algoritmos = {}
    algoritmos['GCG'] = int(input("GCG (0/1): "))
    algoritmos['GRASP'] = int(input("GRASP (0/1): "))
    algoritmos['GRASP_GCG'] = int(input("GRASP_GCG (0/1): "))
    algoritmos['Random_Key'] = int(input("Random Key (0/1): "))
    # Novos algoritmos
    algoritmos['BRKGA'] = int(input("BRKGA (0/1): "))
    algoritmos['GreedyInisol'] = int(input("GreedyInisol (0/1): "))
    algoritmos['GreedySearch'] = int(input("GreedySearch (0/1): "))
    
    num_execucoes = 1
    # Considera algoritmos iterativos
    if (algoritmos['GRASP'] or algoritmos['GRASP_GCG'] or algoritmos['Random_Key'] or algoritmos['BRKGA']):
        while True:
            try:
                num_execucoes = int(input("\nNúmero de execuções para algoritmos iterativos: "))
                if num_execucoes > 0:
                    break
                print("O número de execuções deve ser maior que 0.")
            except ValueError:
                print("Por favor, digite um número válido.")
    
    return {
        'dimensoes': dimensoes_selecionadas,
        'algoritmos': algoritmos,
        'num_execucoes': num_execucoes
    }

def calcular_dimensoes_area(poligonos, dimensoes_percentuais, proporcao_base_altura=1.5):
    """
    Para cada valor percentual informado, calcula as dimensões necessárias para que a área total dos polígonos
    seja igual à área desejada. Ajusta se as dimensões não comportam a maior peça.
    Retorna uma lista com itens do formato: [base, altura, percentual_corrigido]
    """
    area_total = sum(calcular_area(poligono) for poligono in poligonos)
    
    max_dx = 0
    max_dy = 0
    for poligono in poligonos:
        xs = [p[0] for p in poligono]
        ys = [p[1] for p in poligono]
        dx = max(xs) - min(xs)
        dy = max(ys) - min(ys)
        if dx > max_dx:
            max_dx = dx
        if dy > max_dy:
            max_dy = dy

    dimensoes_calculadas = []
    
    def calcular_dimensoes(area):
        altura = (area / proporcao_base_altura) ** 0.5
        base = altura * proporcao_base_altura
        return round(base, 2), round(altura, 2)
    
    for perc in dimensoes_percentuais:
        area_desejada = area_total / (perc / 100.0)
        base, altura = calcular_dimensoes(area_desejada)
        
        base_original, altura_original = base, altura
        perc_original = perc
        
        if base < max_dx or altura < max_dy:
            fator = max(max_dx / base if base > 0 else 1, max_dy / altura if altura > 0 else 1)
            base *= fator
            altura *= fator
            base += 5
            altura += 5
            base = int(base)
            altura = int(altura)
            nova_area = base * altura
            perc = int(area_total / nova_area * 100.0)
            print(f"Aviso: A dimensão calculada para {perc_original}% (Base={base_original}, Altura={altura_original})")
            print(f"não comporta a maior peça (max_dx={max_dx}, max_dy={max_dy}).")
            print(f"Nova dimensão ajustada: Base={base}, Altura={altura}, com percentual corrigido de {perc}%.\n")
        
        dimensoes_calculadas.append([int(base), int(altura), perc])
    
    return dimensoes_calculadas

def save_algorithm_results(algorithm_name, results, tempos, areas, pecas, fechos, boxes, solucoes, num_execucoes, results_file, env):
    with open(results_file, "a", encoding='utf-8') as file:
        file.write("-"*190 + "\n")
        file.write(f"RESULTADOS {algorithm_name}:\n")
        
        for i in range(num_execucoes):
            file.write(results[i] + "\n")
            file.write(f"Peças {algorithm_name} {i+1}: {solucoes[i]}\n\n")
        
        if num_execucoes > 1:
            media_tempo = round(sum(tempos) / num_execucoes, 2)
            desvio_tempo = round((sum((x - media_tempo)**2 for x in tempos) / num_execucoes)**0.5, 2)
            
            media_area = round(sum(areas) / num_execucoes, 2)
            desvio_area = round((sum((x - media_area)**2 for x in areas) / num_execucoes)**0.5, 2)
            
            media_pecas = round(sum(pecas) / num_execucoes, 2)
            desvio_pecas = round((sum((x - media_pecas)**2 for x in pecas) / num_execucoes)**0.5, 2)
            
            media_fecho = round(sum(fechos) / num_execucoes, 2)
            desvio_fecho = round((sum((x - media_fecho)**2 for x in fechos) / num_execucoes)**0.5, 2)
            
            media_box = round(sum(boxes) / num_execucoes, 2)
            desvio_box = round((sum((x - media_box)**2 for x in boxes) / num_execucoes)**0.5, 2)
            
            estatisticas = (
                f"\nEstatísticas {algorithm_name}:\n"
                f"Tempo - Média: {media_tempo}s, Desvio: {desvio_tempo}s\n"
                f"Área - Média: {media_area}%, Desvio: {desvio_area}%\n"
                f"Peças - Média: {media_pecas}, Desvio: {desvio_pecas}\n"
                f"Fecho - Média: {media_fecho}%, Desvio: {desvio_fecho}%\n"
                f"Box - Média: {media_box}%, Desvio: {desvio_box}%"
            )
            file.write(estatisticas + "\n")
        
        file.write("-"*190 + "\n\n")

def save_cycle_results(algorithm_name, dataset, cycle_num, tempo, area, pecas, fecho, box, folder, pecas_posicionadas):
    today = date.today().strftime("%Y-%m-%d")
    txt_file = os.path.join(folder, f"resultados_ciclos_{today}.txt")
    csv_file = os.path.join(folder, f"resultados_ciclos_{today}.csv")
    
    with open(txt_file, "a", encoding='utf-8') as file:
        if cycle_num == 1:
            file.write("-"*100 + "\n")
            file.write(f"Dataset: {dataset} | Algoritmo: {algorithm_name}\n")
            file.write("{:<10} {:<10} {:<10} {:<10} {:<15} {:<15}\n".format(
                "Ciclo", "Tempo(s)", "Área(%)", "Peças", "Fecho(%)", "Box(%)"))
        
        file.write("{:<10} {:<10.2f} {:<10.2f} {:<10} {:<15.2f} {:<15.2f}".format(
            cycle_num, tempo, area, pecas, fecho, box))
        file.write(f" Solução: {pecas_posicionadas}\n")
    
    csv_header = ["Dataset", "Algoritmo", "Ciclo", "Tempo(s)", "Área(%)", "Peças", "Fecho(%)", "Box(%)"]
    file_exists = os.path.isfile(csv_file)
    
    with open(csv_file, "a", newline='', encoding='utf-8') as file:
        writer = csv.DictWriter(file, fieldnames=csv_header)
        if not file_exists:
            writer.writeheader()
        writer.writerow({
            "Dataset": dataset,
            "Algoritmo": algorithm_name,
            "Ciclo": cycle_num,
            "Tempo(s)": f"{tempo:.2f}",
            "Área(%)": f"{area:.2f}",
            "Peças": pecas,
            "Fecho(%)": f"{fecho:.2f}",
            "Box(%)": f"{box:.2f}"
        })

def save_dimension_summary(dataset, dimension, results, filename):
    metrics_mapping = {"Tempo": "tempos", "Área": "areas", "Peças": "pecas", "Fecho": "fechos", "Box": "boxes"}
    
    headers = ["Dataset", "Dimensão"]
    for algo in results.keys():
        for metric in metrics_mapping.keys():
            headers.extend([f"{algo}_{metric}_Media", f"{algo}_{metric}_Desvio", f"{algo}_{metric}_Max"])
    
    row = {"Dataset": dataset, "Dimensão": f"{dimension[0]}x{dimension[1]}"}
    
    for algo, data in results.items():
        for metric, key in metrics_mapping.items():
            values = data.get(key, [])
            media = sum(values)/len(values) if values else 0
            desvio = (sum((x - media)**2 for x in values)/len(values))**0.5 if len(values) > 1 else 0
            max_val = max(values) if values else 0
            
            row[f"{algo}_{metric}_Media"] = round(media, 2)
            row[f"{algo}_{metric}_Desvio"] = round(desvio, 2)
            row[f"{algo}_{metric}_Max"] = round(max_val, 2)
    
    file_exists = os.path.isfile(filename)
    with open(filename, "a", newline='', encoding='utf-8') as file:
        writer = csv.DictWriter(file, fieldnames=headers)
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)

def save_um_ciclo_summary(dataset, dimension, results, filename):
    headers = ["Dataset", "Dimensão"]
    for algo in results.keys():
        headers.extend([f"{algo}_pecas", f"{algo}_%preenchimento", f"{algo}_tempo"])
    
    row = {"Dataset": dataset, "Dimensão": f"{dimension[0]}x{dimension[1]}"}
    for algo, data in results.items():
        pecas_avg = round(sum(data['pecas'])/len(data['pecas']), 2) if data.get('pecas') else 0
        areas_avg = round(sum(data['areas'])/len(data['areas']), 2) if data.get('areas') else 0
        tempos_avg = round(sum(data['tempos'])/len(data['tempos']), 2) if data.get('tempos') else 0
        row[f"{algo}_pecas"] = pecas_avg
        row[f"{algo}_%preenchimento"] = areas_avg
        row[f"{algo}_tempo"] = tempos_avg

    file_exists = os.path.isfile(filename)
    with open(filename, "a", newline='', encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=headers)
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)

def juntar_datasets(datasets):
    name = ''.join(datasets)
    output = []
    lista_final = []
    for dataset in datasets:
        
       
        output.append(ler_poligonos('/home/fsilvestre/Cutting_Stock_Problem/Datasets/Embraer/' + dataset +'.dat'))

    for lista in output:
        for pol in lista:
            lista_final.append(pol)

    tamanho = len(lista_final)

    escrever_poligonos(lista_final, '/home/fsilvestre/Cutting_Stock_Problem/Datasets/Embraer/' + name + '.dat')

    return name




if __name__ == "__main__":
    # Cria a pasta de teste e obtém seu nome
    pasta_teste = criar_pasta_teste()
    
    config = configurar_teste()
    
    results_file = create_results_file(pasta_teste)
    margem = int(input("\nCom margem?(0/5) "))
    uc = int(input("\nUm ciclo?(0/1) "))
    um_ciclo = uc == 1
    mc = int(input("\nMulti-ciclo?(0/1) "))
    multi_ciclo = mc == 1

    separados = int(input("\nDatasets separados ou juntos? (0/1)"))
    separados = separados == 0
    fr = 1

    current_time = datetime.now()

    instancias = ['embraer_55', 'embraer_610', 'embraer_2338', 'embraer_3274',
                      'embraer_3089', 'embraer_3086', 'embraer_3153']
    if not separados:
        datasets = []
        for ins in instancias:
            tag = int(input(f"Adicionar {ins}? (0/1)\n"))
            if tag == 1:
                datasets.append(ins)

        name = juntar_datasets(datasets)
        instancias = [name]
        
        print(instancias)

    
    rotacoes = [0, 1, 2, 3]
    vizinhanca = 25
    q1 = 50
    q2 = 20
    q3 = 30
    pecas_inisol = 0
    max_iter = 20

    pop = 40
    gen = 20

    um_ciclo_csv = os.path.join(pasta_teste, f"um_ciclo_summary_{date.today().strftime('%Y-%m-%d')}.csv")

    for instancia in reversed(instancias):
        print('\n')
        env = CSP(instancia, plot=False, render=False)
        dimensoes = calcular_dimensoes_area(env.nova_lista, config['dimensoes'])
        
        with open(results_file, "a", encoding='utf-8') as file:
            file.write("\n" + "|"*190 + "\n")
            file.write(f"Instância: {instancia}\n")
            file.write(f"Data e Hora: {current_time.strftime('%d/%m/%Y %H:%M')}\n")
            file.write(f"Margem: {margem > 0}\n")
            file.write(str(config) + '\n')
            file.write("|"*190 + "\n\n")

        Stime = time.time()
        tabela_nfps = pre_processar_NFP(rotacoes, env.nova_lista, margem)
        Etime = time.time()
        tempo_nfp = round(Etime - Stime, 2)
        print(f'\nNfps pré-processados em {tempo_nfp}s')

        for idx, (perc, dimensao) in enumerate(zip(config['dimensoes'], dimensoes)):
            print(f"\nTestando para {perc}% -> Dimensão (Base x Altura): {dimensao}")
            env = CSP(instancia, plot=False, render=False, Base=dimensao[0], Altura=dimensao[1])
            env_ciclo = CSP(instancia, plot=False, render=False, Base=int(dimensao[0]*fr), Altura=int(dimensao[1]*fr))
            
            dimension_results = {}
            
            with open(results_file, "a", encoding='utf-8') as file:
                file.write("\n" + "="*190 + "\n")
                file.write(f"Dimensão: {dimensao}\n")
                file.write(f"Máximo: {perc}%\n")
                file.write(f"Tempo pré-processamento NFPs: {tempo_nfp}s\n")              
                file.write("="*190 + "\n\n")

            # GCG
            if config['algoritmos']['GCG']:
                if um_ciclo:
                    Stime = time.time()
                    pecas_GCG, area, lista_ciclos = GCG(tabela_nfps, instancia, plot=False, render=False, 
                                           base=dimensao[0], altura=dimensao[1], margem=margem)
                    Etime = time.time()
                    tempo_GCG = round(Etime - Stime, 2)
                    
                    area_percent = round(area * 100, 2)
                    fecho_percent = round(calcular_fecho_area(pecas_GCG)/env.area * 100, 2)
                    box_percent = round(area_fecho_retangular(pecas_GCG)/env.area * 100, 2)
                    
                    resultado_GCG = (
                        f"GCG: Tempo: {tempo_GCG}s, Area: {area_percent}%, "
                        f"Peças: {len(pecas_GCG)}/{env.max_pecas}, "
                        f"Fecho: {fecho_percent}%, Box: {box_percent}%"
                    )
                    
                    save_algorithm_results(
                        'GCG', 
                        [resultado_GCG], 
                        [tempo_GCG], 
                        [area_percent], 
                        [len(pecas_GCG)], 
                        [fecho_percent], 
                        [box_percent], 
                        [pecas_GCG], 
                        1, 
                        results_file, 
                        env
                    )
                    dimension_results['GCG'] = {
                        'tempos': [tempo_GCG],
                        'areas': [area_percent],
                        'pecas': [len(pecas_GCG)],
                        'fechos': [fecho_percent],
                        'boxes': [box_percent]
                    }
                    print(resultado_GCG)
                
                if multi_ciclo:
                    dimensao[0] *= fr
                    dimensao[1] *= fr
                    arquivo_ciclo = os.path.join(pasta_teste, f"{instancia}_novo_ciclo_gcg.dat")
                    with open(arquivo_ciclo, 'w') as arquivo:
                        pass
                    
                    ciclo_count = 1
                    lista_ciclos = copy.deepcopy(env.nova_lista)
                    
                    Stime = time.time()
                    pecas_GCG, area, lista_ciclos = GCG(tabela_nfps, instancia, plot=False, render=False, 
                                                        base=int(dimensao[0]), altura=int(dimensao[1]), margem=margem)
                    Etime = time.time()
                    tempo_GCG = round(Etime - Stime, 2)
                    area_percent = round(area * 100, 2)
                    fecho_percent = round(calcular_fecho_area(pecas_GCG)/env.area * 100, 2)
                    box_percent = round(area_fecho_retangular(pecas_GCG)/env.area * 100, 2)
                    
                    save_cycle_results(
                        'GCG',
                        instancia,
                        ciclo_count,
                        tempo_GCG,
                        area_percent,
                        len(pecas_GCG),
                        fecho_percent,
                        box_percent,
                        pasta_teste,
                        pecas_GCG
                    )
                    
                    while len(lista_ciclos) > 0:
                        ciclo_count += 1
                        escrever_poligonos(lista_ciclos, arquivo_ciclo)
                        Stime = time.time()
                        pecas_GCG, area, lista_ciclos = GCG(tabela_nfps, arquivo_ciclo, plot=False, render=False, 
                                                            base=int(dimensao[0]), altura=int(dimensao[1]), margem=margem, suavizar=False)
                        Etime = time.time()
                        tempo_GCG = round(Etime - Stime, 2)
                        area_percent = round(area * 100, 2)
                        fecho_percent = round(calcular_fecho_area(pecas_GCG)/env.area * 100, 2)
                        box_percent = round(area_fecho_retangular(pecas_GCG)/env.area * 100, 2)
                        
                        save_cycle_results(
                            'GCG',
                            instancia,
                            ciclo_count,
                            tempo_GCG,
                            area_percent,
                            len(pecas_GCG),
                            fecho_percent,
                            box_percent,
                            pasta_teste,
                            pecas_GCG
                        )
                    dimensao[0] /= fr
                    dimensao[1] /= fr
            
            # GRASP
            if config['algoritmos']['GRASP']:
                if um_ciclo:
                    resultados_GRASP = []
                    tempos_GRASP = []
                    areas_GRASP = []
                    pecas_GRASP = []
                    fechos_GRASP = []
                    boxes_GRASP = []
                    solucoes_GRASP = []

                    for i in range(config['num_execucoes']):
                        Stime = time.time()
                        I1, pecas_G, lista = GRASP(max_iter, vizinhanca, q1, q2, q3, rotacoes,
                                                   pecas_inisol, instancia, tabela_nfps, 
                                                   base=dimensao[0], altura=dimensao[1])
                        Etime = time.time()
                        tempo_G = round(Etime - Stime, 2)
                        
                        area_calc = round((sum(calcular_area(p) for p in pecas_G)/env.area)*100, 2)
                        fecho_calc = round(calcular_fecho_area(pecas_G)/env.area*100, 2)
                        box_calc = round(area_fecho_retangular(pecas_G)/env.area*100, 2)
                        
                        resultado = (
                            f"GRASP: Tempo: {tempo_G}s, Area: {area_calc}%, "
                            f"Peças: {len(pecas_G)}/{env.max_pecas}, "
                            f"Fecho: {fecho_calc}%, Box: {box_calc}%"
                        )
                        
                        resultados_GRASP.append(resultado)
                        tempos_GRASP.append(tempo_G)
                        areas_GRASP.append(area_calc)
                        pecas_GRASP.append(len(pecas_G))
                        fechos_GRASP.append(fecho_calc)
                        boxes_GRASP.append(box_calc)
                        solucoes_GRASP.append(pecas_G)
                        print(resultado)

                    save_algorithm_results(
                        'GRASP', 
                        resultados_GRASP, 
                        tempos_GRASP, 
                        areas_GRASP, 
                        pecas_GRASP, 
                        fechos_GRASP, 
                        boxes_GRASP, 
                        solucoes_GRASP, 
                        config['num_execucoes'], 
                        results_file, 
                        env
                    )
                    dimension_results['GRASP'] = {
                        'tempos': tempos_GRASP,
                        'areas': areas_GRASP,
                        'pecas': pecas_GRASP,
                        'fechos': fechos_GRASP,
                        'boxes': boxes_GRASP
                    }
                
                if multi_ciclo:
                    dimensao[0] *= fr
                    dimensao[1] *= fr
                    instancia_ciclo = f"{instancia}_novo_ciclo_g"
                    arquivo_ciclo = os.path.join(pasta_teste, f"{instancia}_novo_ciclo_g.dat")
                    with open(arquivo_ciclo, 'w') as arquivo:
                        pass
                    
                    ciclo_count = 1

                    Stime = time.time()
                    I1, pecas_G, lista_ciclos = GRASP(max_iter, vizinhanca, q1, q2, q3, rotacoes,
                                                        pecas_inisol, instancia, tabela_nfps, 
                                                        base=int(dimensao[0]), altura=int(dimensao[1]))
                    Etime = time.time()
                    tempo_G = round(Etime - Stime, 2)
                    
                    area_calc = round((sum(calcular_area(p) for p in pecas_G)/env_ciclo.area)*100, 2)
                    fecho_calc = round(calcular_fecho_area(pecas_G)/env_ciclo.area*100, 2)
                    box_calc = round(area_fecho_retangular(pecas_G)/env_ciclo.area*100, 2)
                    
                    save_cycle_results(
                        'GRASP',
                        instancia,
                        ciclo_count,
                        tempo_G,
                        area_calc,
                        len(pecas_G),
                        fecho_calc,
                        box_calc,
                        pasta_teste,
                        pecas_G
                    )
                                            
                    
                    while len(lista_ciclos) > 0:
                        ciclo_count += 1
                        escrever_poligonos(lista_ciclos, arquivo_ciclo)

                        Stime = time.time()
                        I1, pecas_G, lista_ciclos = GRASP(max_iter, vizinhanca, q1, q2, q3, rotacoes,
                                                          pecas_inisol, arquivo_ciclo, tabela_nfps, 
                                                          base=int(dimensao[0]), altura=int(dimensao[1]),suavizar=False)
                        Etime = time.time()
                        tempo_G = round(Etime - Stime, 2)
                        
                        area_calc = round((sum(calcular_area(p) for p in pecas_G)/env_ciclo.area)*100, 2)
                        fecho_calc = round(calcular_fecho_area(pecas_G)/env_ciclo.area*100, 2)
                        box_calc = round(area_fecho_retangular(pecas_G)/env_ciclo.area*100, 2)
                        
                        save_cycle_results(
                            'GRASP',
                            instancia,
                            ciclo_count,
                            tempo_G,
                            area_calc,
                            len(pecas_G),
                            fecho_calc,
                            box_calc,
                            pasta_teste,
                            pecas_G
                        )
                        
                                                  
                    dimensao[0] /= fr
                    dimensao[1] /= fr
            
            # GRASP_GCG
            if config['algoritmos']['GRASP_GCG']:
                if um_ciclo:
                    resultados_GRASP_GCG = []
                    tempos_GRASP_GCG = []
                    areas_GRASP_GCG = []
                    pecas_GRASP_GCG = []
                    fechos_GRASP_GCG = []
                    boxes_GRASP_GCG = []
                    solucoes_GRASP_GCG = []

                    for i in range(config['num_execucoes']):
                        Stime = time.time()
                        I1, pecas_GG,_ = GRASP_GCG(max_iter, vizinhanca, q1, q2, q3, rotacoes,
                                                pecas_inisol, instancia, tabela_nfps, 
                                                base=dimensao[0], altura=dimensao[1])
                        Etime = time.time()
                        tempo_GG = round(Etime - Stime, 2)
                        
                        area_calc = round((sum(calcular_area(p) for p in pecas_GG)/env.area)*100, 2)
                        fecho_calc = round(calcular_fecho_area(pecas_GG)/env.area*100, 2)
                        box_calc = round(area_fecho_retangular(pecas_GG)/env.area*100, 2)
                        
                        resultado = (
                            f"GRASP_GCG: Tempo: {tempo_GG}s, Area: {area_calc}%, "
                            f"Peças: {len(pecas_GG)}/{env.max_pecas}, "
                            f"Fecho: {fecho_calc}%, Box: {box_calc}%"
                        )
                        
                        resultados_GRASP_GCG.append(resultado)
                        tempos_GRASP_GCG.append(tempo_GG)
                        areas_GRASP_GCG.append(area_calc)
                        pecas_GRASP_GCG.append(len(pecas_GG))
                        fechos_GRASP_GCG.append(fecho_calc)
                        boxes_GRASP_GCG.append(box_calc)
                        solucoes_GRASP_GCG.append(pecas_GG)
                        print(resultado)

                    save_algorithm_results(
                        'GRASP_GCG', 
                        resultados_GRASP_GCG, 
                        tempos_GRASP_GCG, 
                        areas_GRASP_GCG, 
                        pecas_GRASP_GCG, 
                        fechos_GRASP_GCG, 
                        boxes_GRASP_GCG, 
                        solucoes_GRASP_GCG, 
                        config['num_execucoes'], 
                        results_file, 
                        env
                    )
                    dimension_results['GRASP_GCG'] = {
                        'tempos': tempos_GRASP_GCG,
                        'areas': areas_GRASP_GCG,
                        'pecas': pecas_GRASP_GCG,
                        'fechos': fechos_GRASP_GCG,
                        'boxes': boxes_GRASP_GCG
                    }
                
                if multi_ciclo:
                    dimensao[0] *= fr
                    dimensao[1] *= fr
                    arquivo_ciclo = os.path.join(pasta_teste, f"{instancia}_novo_ciclo_gg.dat")
                    with open(arquivo_ciclo, 'w') as arquivo:
                        pass
                    
                    ciclo_count = 1
                    lista_ciclos_gg = copy.deepcopy(env.nova_lista)
                    
                    Stime = time.time()
                    I1, pecas_GG, lista_ciclos_gg = GRASP_GCG(max_iter, vizinhanca, q1, q2, q3, rotacoes,
                                                            pecas_inisol, instancia, tabela_nfps, 
                                                            base=int(dimensao[0]), altura=int(dimensao[1]))
                    Etime = time.time()
                    tempo_GG = round(Etime - Stime, 2)
                    area_calc = round((sum(calcular_area(p) for p in pecas_GG)/env.area)*100, 2)
                    fecho_calc = round(calcular_fecho_area(pecas_GG)/env.area*100, 2)
                    box_calc = round(area_fecho_retangular(pecas_GG)/env.area*100, 2)
                    
                    save_cycle_results(
                        'GRASP_GCG',
                        instancia,
                        ciclo_count,
                        tempo_GG,
                        area_calc,
                        len(pecas_GG),
                        fecho_calc,
                        box_calc,
                        pasta_teste,
                        pecas_GG
                    )
                    
                    while len(lista_ciclos_gg) > 0:
                        ciclo_count += 1
                        escrever_poligonos(lista_ciclos_gg, arquivo_ciclo)
                        Stime = time.time()
                        I1, pecas_GG, lista_ciclos_gg = GRASP_GCG(max_iter, vizinhanca, q1, q2, q3, rotacoes,
                                                                pecas_inisol, arquivo_ciclo, tabela_nfps, 
                                                                base=int(dimensao[0]), altura=int(dimensao[1]),suavizar=False)
                        Etime = time.time()
                        tempo_GG = round(Etime - Stime, 2)
                        area_calc = round((sum(calcular_area(p) for p in pecas_GG)/env.area)*100, 2)
                        fecho_calc = round(calcular_fecho_area(pecas_GG)/env.area*100, 2)
                        box_calc = round(area_fecho_retangular(pecas_GG)/env.area*100, 2)
                        
                        save_cycle_results(
                            'GRASP_GCG',
                            instancia,
                            ciclo_count,
                            tempo_GG,
                            area_calc,
                            len(pecas_GG),
                            fecho_calc,
                            box_calc,
                            pasta_teste,
                            pecas_GG
                        )
                    dimensao[0] /= fr
                    dimensao[1] /= fr
            
            # Random Key
            if config['algoritmos']['Random_Key']:
                if um_ciclo:
                    resultados_RK = []
                    tempos_RK = []
                    areas_RK = []
                    pecas_RK = []
                    fechos_RK = []
                    boxes_RK = []
                    solucoes_RK = []

                    for i in range(config['num_execucoes']):
                        Stime = time.time()
                        solucao, lista_ciclos = RKGA(ambiente=env, rotacoes=rotacoes,
                                                     pop_size=pop, n_generations=gen,
                                                     tabela_nfps=tabela_nfps,
                                                     base=dimensao[0], altura=dimensao[1])
                        Etime = time.time()
                        tempo_RK_i = round(Etime - Stime, 2)
                        
                        ambiente = CSP(dataset=instancia, render=False, plot=False,
                                       Base=dimensao[0], Altura=dimensao[1])
                        x, y = ambiente.cordenadas_area[3]
                        index = ambiente.lista.index(solucao[0][-1])
                        ambiente.acao(index, x, y, 0, False)
                        
                        for peca in solucao[1:]:
                            if peca[-1] in ambiente.lista:
                                ambiente.acao(peca[0], peca[1], peca[2], peca[3], peca[4], True)
                        
                        area_calc = round((ambiente.area_ocupada/ambiente.area)*100, 2)
                        fecho_calc = round(calcular_fecho_area(ambiente.pecas_posicionadas)/env.area*100, 2)
                        box_calc = round(area_fecho_retangular(ambiente.pecas_posicionadas)/env.area*100, 2)
                        
                        resultado = (
                            f"Random Key Execução {i+1}: Tempo: {tempo_RK_i}s, Area: {area_calc}%, "
                            f"Peças: {len(ambiente.pecas_posicionadas)}/{env.max_pecas}, "
                            f"Fecho: {fecho_calc}%, Box: {box_calc}%"
                        )
                        
                        resultados_RK.append(resultado)
                        tempos_RK.append(tempo_RK_i)
                        areas_RK.append(area_calc)
                        pecas_RK.append(len(ambiente.pecas_posicionadas))
                        fechos_RK.append(fecho_calc)
                        boxes_RK.append(box_calc)
                        solucoes_RK.append(ambiente.pecas_posicionadas)
                        print(resultado)

                    save_algorithm_results(
                        'Random Key', 
                        resultados_RK, 
                        tempos_RK, 
                        areas_RK, 
                        pecas_RK, 
                        fechos_RK, 
                        boxes_RK, 
                        solucoes_RK, 
                        config['num_execucoes'], 
                        results_file, 
                        env
                    )
                    dimension_results['Random_Key'] = {
                        'tempos': tempos_RK,
                        'areas': areas_RK,
                        'pecas': pecas_RK,
                        'fechos': fechos_RK,
                        'boxes': boxes_RK
                    }
                if multi_ciclo:
                    dimensao[0] *= fr
                    dimensao[1] *= fr
                    arquivo_ciclo = os.path.join(pasta_teste, f"{instancia}_novo_ciclo_rk.dat")
                    with open(arquivo_ciclo, 'w') as arquivo:
                        pass
                    
                    ciclo_count = 1
                    lista_ciclos_rk = copy.deepcopy(env.nova_lista)
                    
                    Stime = time.time()
                    solucao, lista_ciclos_rk = RKGA(ambiente=env, rotacoes=rotacoes,
                                                                pop_size=pop, n_generations=gen,
                                                                tabela_nfps=tabela_nfps,
                                                                base=int(dimensao[0]), altura=int(dimensao[1]))
                    Etime = time.time()
                    tempo_RK = round(Etime - Stime, 2)
                    ambiente = CSP(dataset=instancia, render=False, plot=False,
                                    Base=int(dimensao[0]), Altura=int(dimensao[1]))
                    x, y = ambiente.cordenadas_area[3]
                    index = ambiente.lista.index(solucao[0][-1])
                    ambiente.acao(index, x, y, 0, False)
                    
                    for peca in solucao[1:]:
                        if peca[-1] in ambiente.lista:
                            ambiente.acao(peca[0], peca[1], peca[2], peca[3], peca[4], True)
                    
                    area_calc = round((ambiente.area_ocupada/ambiente.area)*100, 2)
                    fecho_calc = round(calcular_fecho_area(ambiente.pecas_posicionadas)/env.area*100, 2)
                    box_calc = round(area_fecho_retangular(ambiente.pecas_posicionadas)/env.area*100, 2)
                    
                    save_cycle_results(
                        'Random Key',
                        instancia,
                        ciclo_count,
                        tempo_RK,
                        area_calc,
                        len(ambiente.pecas_posicionadas),
                        fecho_calc,
                        box_calc,
                        pasta_teste,
                        ambiente.pecas_posicionadas
                    )
                    
                    while len(lista_ciclos_rk) > 0:
                        
                        ciclo_count += 1
                        escrever_poligonos(lista_ciclos_rk, arquivo_ciclo)
                        env_temp = CSP(arquivo_ciclo, plot=False, render=False, Base=int(dimensao[0]*fr), Altura=int(dimensao[1]*fr),suavizar=False)

                        Stime = time.time()
                        solucao, lista_ciclos_rk = RKGA(ambiente=env_temp, rotacoes=rotacoes,
                                                                    pop_size=pop, n_generations=gen,
                                                                    tabela_nfps=tabela_nfps,
                                                                    base=int(dimensao[0]), altura=int(dimensao[1]), suavizar=False)
                        Etime = time.time()
                        tempo_RK = round(Etime - Stime, 2)
                        ambiente = CSP(dataset=arquivo_ciclo, render=False, plot=False,
                                        Base=dimensao[0], Altura=dimensao[1],suavizar=False)
                        x, y = ambiente.cordenadas_area[3]
                        index = ambiente.lista.index(solucao[0][-1])
                        ambiente.acao(index, x, y, 0, False)
                        
                        for peca in solucao[1:]:
                            if peca[-1] in ambiente.lista:
                                ambiente.acao(peca[0], peca[1], peca[2], peca[3], peca[4], True)
                        
                        area_calc = round((ambiente.area_ocupada/ambiente.area)*100, 2)
                        fecho_calc = round(calcular_fecho_area(ambiente.pecas_posicionadas)/env.area*100, 2)
                        box_calc = round(area_fecho_retangular(ambiente.pecas_posicionadas)/env.area*100, 2)
                        
                        save_cycle_results(
                            'Random Key',
                            instancia,
                            ciclo_count,
                            tempo_RK,
                            area_calc,
                            len(ambiente.pecas_posicionadas),
                            fecho_calc,
                            box_calc,
                            pasta_teste,
                            ambiente.pecas_posicionadas
                        )
                    dimensao[0] /= fr
                    dimensao[1] /= fr

            # BRKGA (saída igual ao RKGA)
            if config['algoritmos']['BRKGA']:
                if um_ciclo:
                    resultados_BRKGA = []
                    tempos_BRKGA = []
                    areas_BRKGA = []
                    pecas_BRKGA = []
                    fechos_BRKGA = []
                    boxes_BRKGA = []
                    solucoes_BRKGA = []
                    for i in range(config['num_execucoes']):
                        Stime = time.time()
                        solucao, lista_ciclos = BRKGA(ambiente=env, rotacoes=rotacoes,
                                        pop_size=pop, n_generations=gen,
                                        tabela_nfps=tabela_nfps,
                                        base=dimensao[0], altura=dimensao[1])
                        Etime = time.time()
                        tempo_BRKGA = round(Etime - Stime, 2)
                        ambiente = CSP(dataset=instancia, render=False, plot=False,
                                       Base=dimensao[0], Altura=dimensao[1])
                        x, y = ambiente.cordenadas_area[3]
                        index = ambiente.lista.index(solucao[0][-1])
                        ambiente.acao(index, x, y, 0, False)
                        for peca in solucao[1:]:
                            if peca[-1] in ambiente.lista:
                                ambiente.acao(peca[0], peca[1], peca[2], peca[3], peca[4], True)
                        area_calc = round((ambiente.area_ocupada/ambiente.area)*100, 2)
                        fecho_calc = round(calcular_fecho_area(ambiente.pecas_posicionadas)/env.area*100, 2)
                        box_calc = round(area_fecho_retangular(ambiente.pecas_posicionadas)/env.area*100, 2)
                        resultado = (
                            f"BRKGA Execução {i+1}: Tempo: {tempo_BRKGA}s, Area: {area_calc}%, "
                            f"Peças: {len(ambiente.pecas_posicionadas)}/{env.max_pecas}, "
                            f"Fecho: {fecho_calc}%, Box: {box_calc}%"
                        )
                        resultados_BRKGA.append(resultado)
                        tempos_BRKGA.append(tempo_BRKGA)
                        areas_BRKGA.append(area_calc)
                        pecas_BRKGA.append(len(ambiente.pecas_posicionadas))
                        fechos_BRKGA.append(fecho_calc)
                        boxes_BRKGA.append(box_calc)
                        solucoes_BRKGA.append(ambiente.pecas_posicionadas)
                        print(resultado)
                    
                    save_algorithm_results(
                        'BRKGA',
                        resultados_BRKGA,
                        tempos_BRKGA,
                        areas_BRKGA,
                        pecas_BRKGA,
                        fechos_BRKGA,
                        boxes_BRKGA,
                        solucoes_BRKGA,
                        config['num_execucoes'],
                        results_file,
                        env
                    )
                    dimension_results['BRKGA'] = {
                        'tempos': tempos_BRKGA,
                        'areas': areas_BRKGA,
                        'pecas': pecas_BRKGA,
                        'fechos': fechos_BRKGA,
                        'boxes': boxes_BRKGA
                    }
                if multi_ciclo:
                    dimensao[0] *= fr
                    dimensao[1] *= fr
                    arquivo_ciclo = os.path.join(pasta_teste, f"{instancia}_novo_ciclo_brkga.dat")
                    with open(arquivo_ciclo, 'w') as arquivo:
                        pass
                    ciclo_count = 1
                    lista_ciclos_brkga = copy.deepcopy(env.nova_lista)
                    Stime = time.time()
                    solucao, lista_ciclos_brkga = BRKGA(ambiente=env, rotacoes=rotacoes,
                                                        pop_size=pop, n_generations=gen,
                                                        tabela_nfps=tabela_nfps,
                                                        base=int(dimensao[0]), altura=int(dimensao[1]))
                    Etime = time.time()
                    tempo_BRKGA = round(Etime - Stime, 2)
                    ambiente = CSP(dataset=instancia, render=False, plot=False,
                                   Base=int(dimensao[0]), Altura=int(dimensao[1]))
                    x, y = ambiente.cordenadas_area[3]
                    index = ambiente.lista.index(solucao[0][-1])
                    ambiente.acao(index, x, y, 0, False)
                    for peca in solucao[1:]:
                        if peca[-1] in ambiente.lista:
                            ambiente.acao(peca[0], peca[1], peca[2], peca[3], peca[4], True)
                    area_calc = round((ambiente.area_ocupada/ambiente.area)*100, 2)
                    fecho_calc = round(calcular_fecho_area(ambiente.pecas_posicionadas)/env.area*100, 2)
                    box_calc = round(area_fecho_retangular(ambiente.pecas_posicionadas)/env.area*100, 2)
                    save_cycle_results(
                        'BRKGA',
                        instancia,
                        ciclo_count,
                        tempo_BRKGA,
                        area_calc,
                        len(ambiente.pecas_posicionadas),
                        fecho_calc,
                        box_calc,
                        pasta_teste,
                        ambiente.pecas_posicionadas
                    )
                    while len(lista_ciclos_brkga) > 0:
                        ciclo_count += 1
                        escrever_poligonos(lista_ciclos_brkga, arquivo_ciclo)
                        env_temp = CSP(arquivo_ciclo, plot=False, render=False, Base=int(dimensao[0]*fr), Altura=int(dimensao[1]*fr),suavizar=False)

                        Stime = time.time()
                        solucao, lista_ciclos_brkga = BRKGA(ambiente=env_temp, rotacoes=rotacoes,
                                                            pop_size=pop, n_generations=gen,
                                                            tabela_nfps=tabela_nfps,
                                                            base=int(dimensao[0]), altura=int(dimensao[1]), suavizar=False)
                        Etime = time.time()
                        tempo_BRKGA = round(Etime - Stime, 2)
                        ambiente = CSP(dataset=arquivo_ciclo, render=False, plot=False,
                                       Base=dimensao[0], Altura=dimensao[1],suavizar=False)
                        x, y = ambiente.cordenadas_area[3]
                        index = ambiente.lista.index(solucao[0][-1])
                        ambiente.acao(index, x, y, 0, False)
                        for peca in solucao[1:]:
                            if peca[-1] in ambiente.lista:
                                ambiente.acao(peca[0], peca[1], peca[2], peca[3], peca[4], True)
                        area_calc = round((ambiente.area_ocupada/ambiente.area)*100, 2)
                        fecho_calc = round(calcular_fecho_area(ambiente.pecas_posicionadas)/env.area*100, 2)
                        box_calc = round(area_fecho_retangular(ambiente.pecas_posicionadas)/env.area*100, 2)
                        save_cycle_results(
                            'BRKGA',
                            instancia,
                            ciclo_count,
                            tempo_BRKGA,
                            area_calc,
                            len(ambiente.pecas_posicionadas),
                            fecho_calc,
                            box_calc,
                            pasta_teste,
                            ambiente.pecas_posicionadas
                        )
                    dimensao[0] /= fr
                    dimensao[1] /= fr
            
            # GreedyInisol (saída igual ao GCG)
            if config['algoritmos']['GreedyInisol']:
                if um_ciclo:
                    Stime = time.time()
                    pecas_GI, area, lista_ciclos = GreedyInisol(tabela_nfps, instancia, plot=False, render=False,
                                                 base=dimensao[0], altura=dimensao[1], margem=margem)
                    Etime = time.time()
                    tempo_GI = round(Etime - Stime, 2)
                    area_percent = round(area * 100, 2)
                    fecho_percent = round(calcular_fecho_area(pecas_GI)/env.area*100, 2)
                    box_percent = round(area_fecho_retangular(pecas_GI)/env.area*100, 2)
                    resultado_GI = (
                        f"GreedyInisol: Tempo: {tempo_GI}s, Area: {area_percent}%, "
                        f"Peças: {len(pecas_GI)}/{env.max_pecas}, "
                        f"Fecho: {fecho_percent}%, Box: {box_percent}%"
                    )
                    save_algorithm_results(
                        'GreedyInisol', 
                        [resultado_GI], 
                        [tempo_GI], 
                        [area_percent], 
                        [len(pecas_GI)], 
                        [fecho_percent], 
                        [box_percent], 
                        [pecas_GI], 
                        1, 
                        results_file, 
                        env
                    )
                    dimension_results['GreedyInisol'] = {
                        'tempos': [tempo_GI],
                        'areas': [area_percent],
                        'pecas': [len(pecas_GI)],
                        'fechos': [fecho_percent],
                        'boxes': [box_percent]
                    }
                    print(resultado_GI)
                if multi_ciclo:
                    dimensao[0] *= fr
                    dimensao[1] *= fr
                    arquivo_ciclo = os.path.join(pasta_teste, f"{instancia}_novo_ciclo_gi.dat")
                    with open(arquivo_ciclo, 'w') as arquivo:
                        pass
                    ciclo_count = 1
                    lista_ciclos_gi = copy.deepcopy(env.nova_lista)
                    Stime = time.time()
                    pecas_GI, area, lista_ciclos_gi = GreedyInisol(tabela_nfps, instancia, plot=False, render=False,
                                                                  base=int(dimensao[0]), altura=int(dimensao[1]), margem=margem)
                    Etime = time.time()
                    tempo_GI = round(Etime - Stime, 2)
                    area_percent = round(area * 100, 2)
                    fecho_percent = round(calcular_fecho_area(pecas_GI)/env.area*100, 2)
                    box_percent = round(area_fecho_retangular(pecas_GI)/env.area*100, 2)
                    save_cycle_results(
                        'GreedyInisol',
                        instancia,
                        ciclo_count,
                        tempo_GI,
                        area_percent,
                        len(pecas_GI),
                        fecho_percent,
                        box_percent,
                        pasta_teste,
                        pecas_GI
                    )
                    while len(lista_ciclos_gi) > 0:
                        ciclo_count += 1
                        escrever_poligonos(lista_ciclos_gi, arquivo_ciclo)
                        Stime = time.time()
                        pecas_GI, area, lista_ciclos_gi = GreedyInisol(tabela_nfps, arquivo_ciclo, plot=False, render=False,
                                                                      base=int(dimensao[0]), altura=int(dimensao[1]), margem=margem, suavizar=False)
                        Etime = time.time()
                        tempo_GI = round(Etime - Stime, 2)
                        area_percent = round(area * 100, 2)
                        fecho_percent = round(calcular_fecho_area(pecas_GI)/env.area*100, 2)
                        box_percent = round(area_fecho_retangular(pecas_GI)/env.area*100, 2)
                        save_cycle_results(
                            'GreedyInisol',
                            instancia,
                            ciclo_count,
                            tempo_GI,
                            area_percent,
                            len(pecas_GI),
                            fecho_percent,
                            box_percent,
                            pasta_teste,
                            pecas_GI
                        )
                    dimensao[0] /= fr
                    dimensao[1] /= fr
            
            # GreedySearch (saída igual ao GCG)
            if config['algoritmos']['GreedySearch']:
                if um_ciclo:
                    Stime = time.time()
                    pecas_GS, area, lista_ciclos = GreedySearch(tabela_nfps, instancia, plot=False, render=False,
                                                 base=dimensao[0], altura=dimensao[1], margem=margem)
                    Etime = time.time()
                    tempo_GS = round(Etime - Stime, 2)
                    area_percent = round(area * 100, 2)
                    fecho_percent = round(calcular_fecho_area(pecas_GS)/env.area*100, 2)
                    box_percent = round(area_fecho_retangular(pecas_GS)/env.area*100, 2)
                    resultado_GS = (
                        f"GreedySearch: Tempo: {tempo_GS}s, Area: {area_percent}%, "
                        f"Peças: {len(pecas_GS)}/{env.max_pecas}, "
                        f"Fecho: {fecho_percent}%, Box: {box_percent}%"
                    )
                    save_algorithm_results(
                        'GreedySearch', 
                        [resultado_GS], 
                        [tempo_GS], 
                        [area_percent], 
                        [len(pecas_GS)], 
                        [fecho_percent], 
                        [box_percent], 
                        [pecas_GS], 
                        1, 
                        results_file, 
                        env
                    )
                    dimension_results['GreedySearch'] = {
                        'tempos': [tempo_GS],
                        'areas': [area_percent],
                        'pecas': [len(pecas_GS)],
                        'fechos': [fecho_percent],
                        'boxes': [box_percent]
                    }
                    print(resultado_GS)
                if multi_ciclo:
                    dimensao[0] *= fr
                    dimensao[1] *= fr
                    arquivo_ciclo = os.path.join(pasta_teste, f"{instancia}_novo_ciclo_gs.dat")
                    with open(arquivo_ciclo, 'w') as arquivo:
                        pass
                    ciclo_count = 1
                    lista_ciclos_gs = copy.deepcopy(env.nova_lista)
                    Stime = time.time()
                    pecas_GS, area, lista_ciclos_gs = GreedySearch(tabela_nfps, instancia, plot=False, render=False,
                                                                  base=int(dimensao[0]), altura=int(dimensao[1]), margem=margem)
                    Etime = time.time()
                    tempo_GS = round(Etime - Stime, 2)
                    area_percent = round(area * 100, 2)
                    fecho_percent = round(calcular_fecho_area(pecas_GS)/env.area*100, 2)
                    box_percent = round(area_fecho_retangular(pecas_GS)/env.area*100, 2)
                    save_cycle_results(
                        'GreedySearch',
                        instancia,
                        ciclo_count,
                        tempo_GS,
                        area_percent,
                        len(pecas_GS),
                        fecho_percent,
                        box_percent,
                        pasta_teste,
                        pecas_GS
                    )
                    while len(lista_ciclos_gs) > 0:
                        ciclo_count += 1
                        escrever_poligonos(lista_ciclos_gs, arquivo_ciclo)
                        Stime = time.time()
                        pecas_GS, area, lista_ciclos_gs = GreedySearch(tabela_nfps, arquivo_ciclo, plot=False, render=False,
                                                                      base=int(dimensao[0]), altura=int(dimensao[1]), margem=margem, suavizar=False)
                        Etime = time.time()
                        tempo_GS = round(Etime - Stime, 2)
                        area_percent = round(area * 100, 2)
                        fecho_percent = round(calcular_fecho_area(pecas_GS)/env.area*100, 2)
                        box_percent = round(area_fecho_retangular(pecas_GS)/env.area*100, 2)
                        save_cycle_results(
                            'GreedySearch',
                            instancia,
                            ciclo_count,
                            tempo_GS,
                            area_percent,
                            len(pecas_GS),
                            fecho_percent,
                            box_percent,
                            pasta_teste,
                            pecas_GS
                        )
                    dimensao[0] /= fr
                    dimensao[1] /= fr
            
            # Salvar sumário da dimensão completo (CSV)
            if dimension_results:
                today = date.today().strftime("%Y-%m-%d")
                summary_file = os.path.join(pasta_teste, f"summary_{today}.csv")
                save_dimension_summary(instancia, dimensao, dimension_results, summary_file)
                
                if um_ciclo:
                    save_um_ciclo_summary(instancia, dimensao, dimension_results, um_ciclo_csv)
