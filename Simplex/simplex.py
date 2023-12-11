import easygui
import numpy as np
    
def indices_variaveis_basicas_e_nao_basicas(tableau):
    m, n = tableau.shape

    indices_basicas = []
    indices_nao_basicas = []

    for j in range(n - 1):  # ignora a última coluna (b)
        coluna = tableau[:, j]
        if np.count_nonzero(coluna == 1) == 1 and np.count_nonzero(coluna == 0) == m - 1:
            # encontrou uma coluna com um 1 e o resto zero
            indices_basicas.append(j)
        else:
            # não atendeu às condições de variável básica, adiciona na não básica
            indices_nao_basicas.append(j)

    return np.array(indices_basicas), np.array(indices_nao_basicas)

def novo_cr(c_R, c_B, B_inversa, R): # atualiza as variaveis não basicas na FO
    cr = c_R - c_B @ B_inversa @ R
    return cr

def calcular_b(B_inversa, b): # define o novo valor das variaveis basicas
    b = B_inversa @ b
    return b

def calcular_r(B_inversa, R): # calcula os coeficientes das variaveis não básicas
    r = B_inversa @ R
    return r

def valor_fo(c_B, B_inversa, b):
    fo = c_B @ B_inversa @ b
    return fo

def escolher_variavel_entrada(coeficientes_objetivo):
    # Escolhe a variável que entrará na base
    variavel_entrada = np.argmin(coeficientes_objetivo)
    return variavel_entrada

def substituir_valor(lista, valor_antigo, valor_novo): # trocar valores dentro de uma matriz, função usada para trocar indices das variaveis basicas e nao basicas
    indice = np.where(lista == valor_antigo)
    lista[indice] = valor_novo

def escolher_variavel_saida(matriz_tableau, variavel_entrada):
    m, n = matriz_tableau.shape
    razoes = []

    for i in range(m - 1):  # ignora a ultima linha (fo)
        if matriz_tableau[i, variavel_entrada] > 0:
            razoes.append(matriz_tableau[i, -1] / matriz_tableau[i, variavel_entrada])
        else:
            razoes.append(np.inf)

    variavel_saida = np.argmin(razoes)
    return variavel_saida

def forma_padrao(matriz_tableau, sinal_restricao):
    m, n = matriz_tableau.shape

    # matriz de variáveis de folga ou excesso
    variaveis_adicionais = np.zeros((m, m))

    for i in range(m-1):
        if sinal_restricao[i+1] == '>=':
            variaveis_adicionais[i, i] = -1  # variável de excesso
        elif sinal_restricao[i+1] == '<=':
            variaveis_adicionais[i, i] = 1  # variável de folga

    # remove as colunas que possuem apenas zero, sao colunas de variaveis que n precisaram ser inseridas na matriz
    colunas_nao_zeros = np.any(variaveis_adicionais != 0, axis=0)
    variaveis_adicionais = variaveis_adicionais[:, colunas_nao_zeros]

    # atualiza a matriz_tableau
    matriz_tableau_padrao = np.column_stack((matriz_tableau[:, :-1], variaveis_adicionais, matriz_tableau[:, -1]))

    return matriz_tableau_padrao

def get_input():
    n = int(easygui.enterbox("Insira o numero de variaveis:"))
    m = int(easygui.enterbox("Insira o numero de restricoes:"))

    coeficientes_objetivo = easygui.multenterbox(
        msg="Insira os coeficientes da funcao objetivo:",
        fields=[f"x{i+1}" for i in range(n)]
    )

    choices = ["Max", "Min"]
    msg ="Escolha o tipo do PPL"
    title = "Maximização ou Minimização"
    choice = easygui.choicebox(msg, title, choices)

    if choice == "Max":
        coeficientes_objetivo = [-float(coef) for coef in coeficientes_objetivo]
    else:
        coeficientes_objetivo = [float(coef) for coef in coeficientes_objetivo]

    coeficientes_restricoes = []
    sinal_restricao = {}
    for i in range(m):
        restricao = easygui.multenterbox(
            msg=f"Insira os coeficientes da restricao {i+1}:",
            fields=[f"x{i+1}" for i in range(n)],
        )
        choices = [">=", "<=", "="]
        msg ="Escolha o Sinal da restrição"
        title = "Sinal da restrição"
        choice = easygui.choicebox(msg, title, choices)
        sinal_restricao[i+1] = choice
        restricao = [float(coef) for coef in restricao]
        coeficientes_restricoes.append(restricao)
    
    valor_variaveis_basicas = easygui.multenterbox(
        msg="Insira os valores das restrições:",
        fields=[f"restrição {i+1}" for i in range(m)],
    )
    valor_variaveis_basicas = [[float(valor)] for valor in valor_variaveis_basicas]

    matriz_tableau = np.array(coeficientes_restricoes)
    b = np.array(valor_variaveis_basicas)
    c = np.array(coeficientes_objetivo)
    c = np.append(c, 0)
    matriz_tableau = np.column_stack((matriz_tableau, b))
    matriz_tableau = np.row_stack((matriz_tableau, c))

    return matriz_tableau, b, sinal_restricao

def separar_matrizes_base_nao_base(indices_base, indices_nao_base, matriz_tableau):
    matriz_tableau = np.array(matriz_tableau)

    matriz_base = matriz_tableau[:-1, indices_base]

    matriz_nao_base = matriz_tableau[:-1, indices_nao_base]

    fo = matriz_tableau[-1,:-1]

    c_R = fo[indices_nao_base]
    c_B = fo[indices_base]

    return matriz_base, matriz_nao_base , c_R, c_B

def verifica_solucao_otima(matriz_c):
    return np.any(matriz_c < 0)

def exibir_solucao(indices_base, b, fo):
    mensagem = "Resultado\n"

    for i in range(len(indices_base)):
        mensagem += f"x{indices_base[i]} = {b[i, -1]}\n"
    
    mensagem += f"Valor FO: {fo}"

    easygui.msgbox(mensagem, title="Resultado")

# BIG M
def conta_colunas_com_menos_um(tableau): # função para contar quantas variaveis artificiais serão adicionadas
    linhas, colunas = tableau.shape
    special_columns = 0
    index_linhas = []

    for linha in range(linhas):
        for coluna in range(colunas - 1):  # Excluindo a última coluna (b)
            if tableau[linha, coluna] == -1:
                count_zero = np.count_nonzero(tableau[:, coluna] == 0)
                if count_zero == linhas - 1:
                    special_columns += 1
                    index_linhas.append(linha) # guarda o index da linha que possui o -1 para ser usado mais tarde na função que ira adicionar as colunas das variaveis artificias

    return special_columns, index_linhas

def calcular_m(tableau):
    maior_valor = np.max(np.abs(tableau[-1, :-1]))
    return 10 * maior_valor

def add_variaveis_artificiais(tableau, num_artificial, index_linhas):
    linhas, colunas = tableau.shape
    M = calcular_m(tableau)

    for i in range(num_artificial):
        linha = index_linhas[i]
        # adiciona a nova coluna da variável artificial com 1 na posição correspondente
        tableau = np.insert(tableau, -1, np.zeros(linhas), axis=1)
        tableau[linha, -2] = 1
        tableau[-1, -2] = M

    return tableau, M

def remover_artificial_base(tableau, M): # era pra zerar a linha da f.o das variaveis artificiais pra começar a tira-las da base
    artificial_columns = np.where(tableau[-1, :] == M)[0]

    for coluna in artificial_columns:
        row_index = np.where(tableau[:, coluna] == 1)[0][0]
        tableau[-1, :] -= tableau[row_index, :] * tableau[-1, coluna]

    return tableau

def Big_M(tableau,colunas_com_menos_1,index):
    tableau, M = add_variaveis_artificiais(tableau,colunas_com_menos_1,index)
    tableau = remover_artificial_base(tableau,M)
    return tableau

def exibir_mensagem(mensagem):
    easygui.msgbox(mensagem, title="")

def solve_simplex():
    
    tableau, b, sinal_restricoes = get_input()

    tableau = forma_padrao(tableau,sinal_restricoes)

    # verificando se precisa de variaveis artificiais
    colunas_com_menos_1, index = conta_colunas_com_menos_um(tableau)
    if colunas_com_menos_1 > 0:
        tableau = Big_M(tableau,colunas_com_menos_1,index)

    #selecionar as variaveis que estao atualmente na base
    indices_base, indices_nao_base = indices_variaveis_basicas_e_nao_basicas(tableau)
    
    # escolhe a variavel q irá entrar e retorna seu indice
    variavel_entrada = escolher_variavel_entrada(tableau[-1,indices_nao_base])

    # escolhe a variavel que vai sair da base e retorna seu indice
    variavel_saida = escolher_variavel_saida(tableau,variavel_entrada)

    if variavel_saida == np.inf:
        exibir_mensagem("Solução Ilimitada")
        return 
    
    variavel_saida = indices_base[variavel_saida] # a função escolher_variavel_saida me retorna o indice da variavel que vai sair, porem o indice esta em relação ao vetor das variaveis que estao na base

    # vou na matriz indices_base, encontro o indice do cara que vai sair (variavel_saida) e troco pelo indice da variavel que vai
    # entrar

    substituir_valor(indices_base,variavel_saida,variavel_entrada)

    substituir_valor(indices_nao_base,variavel_entrada,variavel_saida)

    tableau_temp = np.copy(tableau)

    # variavel pra impedir que entre em loop infinito
    # vai guardar o indice da variavel que esta saindo agr, e comparar no proximo loop se é igual a que foi selecionada para entrar
    controle_multiplas = 0 

    while True:
        # divide a matriz A em outras duas matrizes, ja realizando o "pivotamento"
        # pois o recorte sera feito com base no indice passado das variaveis que estao ou não na base
        matriz_B, matriz_R, c_R, c_B = separar_matrizes_base_nao_base(indices_base, indices_nao_base, tableau)

        try:
            B_inversa = np.linalg.inv(matriz_B)
        except np.linalg.LinAlgError:
            exibir_mensagem("A matriz B é singular. O método simplex não pode continuar.")
            return

        c_R = novo_cr(c_R,c_B,B_inversa,matriz_R)

        # atualizar o valor na fo das variaveis nao basicas no tableau temp
        tableau_temp[-1,indices_nao_base] = c_R

        # novo valor das variáveis básicas
        novo_b = calcular_b(B_inversa,b)

        # Calcular novo valor da FO
        fo = valor_fo(c_B,B_inversa,b)
        # atualizar o valor da FO no tableau temp
        tableau_temp[-1,-1] = fo

        # testar se ja estou na solução ótima
        if not verifica_solucao_otima(c_R):
            exibir_solucao(indices_base,novo_b,fo)

            # testar se é multiplas soluções
            if 0 in c_R:
                exibir_mensagem("Multiplas Soluções")
                variavel_entrada_mult = indices_nao_base[np.where(c_R == 0)[0][0]]
                variavel_saida_mult = escolher_variavel_saida(tableau_temp, variavel_entrada_mult)
                variavel_saida_mult = indices_base[variavel_saida_mult]
                controle_multiplas = variavel_saida_mult
                if controle_multiplas == variavel_entrada_mult:
                    break
                substituir_valor(indices_base, variavel_saida_mult, variavel_entrada_mult)
                substituir_valor(indices_nao_base, variavel_entrada_mult, variavel_saida_mult)
                continue
            break

        
        matriz_R = calcular_r(B_inversa,matriz_R)

        # colocando no tableau apenas nas colunas especificas
        tableau_temp[:-1,indices_nao_base] = matriz_R


        variavel_entrada = escolher_variavel_entrada(tableau_temp[-1,indices_nao_base])
        variavel_saida = escolher_variavel_saida(tableau_temp,variavel_entrada)

        if variavel_saida == np.inf:
            exibir_mensagem("Solução Ilimitada")
            return 
        
        variavel_saida = indices_base[variavel_saida]

        # trocando os indices das variaveis basicas e nao basicas
        # para quando o loop reiniciar separar as matrizes com essas referencias
        substituir_valor(indices_base,variavel_saida,variavel_entrada)
        substituir_valor(indices_nao_base,variavel_entrada,variavel_saida)

np.set_printoptions(suppress=True, precision=4)
exibir_mensagem("Big M nao funciona!")
solve_simplex()