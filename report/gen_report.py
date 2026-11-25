# -*- coding: utf-8 -*-
"""Gera o Relatório Final de Engenharia de Software para o projeto Smiles Validator."""
import os
import sys
from docx import Document
from docx.shared import Pt, RGBColor, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH

_HERE = os.path.dirname(os.path.abspath(__file__))
_PIPELINE_PNG = os.path.join(_HERE, "pipeline.png")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(_HERE), "Modelo Estrutural para Relatório Final do Projeto.docx")

doc = Document()

normal = doc.styles["Normal"]
normal.font.name = "Calibri"
normal.font.size = Pt(11)

def h(txt, level=1):
    return doc.add_heading(txt, level=level)

def para(txt="", bold=False, italic=False, size=None):
    p = doc.add_paragraph()
    r = p.add_run(txt)
    r.bold = bold
    r.italic = italic
    if size:
        r.font.size = Pt(size)
    return p

def bullet(txt, bold_prefix=None):
    p = doc.add_paragraph(style="List Bullet")
    if bold_prefix:
        r = p.add_run(bold_prefix); r.bold = True
        p.add_run(txt)
    else:
        p.add_run(txt)
    return p

def numbered(txt):
    p = doc.add_paragraph(style="List Number")
    p.add_run(txt)
    return p

def code_block(lines):
    for ln in lines:
        p = doc.add_paragraph()
        r = p.add_run(ln)
        r.font.name = "Consolas"; r.font.size = Pt(9)
        r.font.color.rgb = RGBColor(0x1F, 0x23, 0x28)
        p.paragraph_format.space_after = Pt(0)
        p.paragraph_format.space_before = Pt(0)

# =========================================================================
# TÍTULO
# =========================================================================
doc.add_heading("Relatório Final do Projeto", level=0)
sub = para("Smiles Validator — Validador Sintático e Químico de Cadeias SMILES", bold=True, size=14)
sub.alignment = WD_ALIGN_PARAGRAPH.CENTER
doc.add_paragraph()

# =========================================================================
# 1. BREVE DESCRIÇÃO
# =========================================================================
h("Breve Descrição", 1)

para("O Smiles Validator é uma biblioteca em Python para validação de cadeias no "
     "formato SMILES (Simplified Molecular Input Line Entry System), uma das notação "
     "textuais mais usada para representar estruturas moleculares. O programa recebe "
     "uma string SMILES e responde se ela é, ao mesmo tempo, sintaticamente bem "
     "formada (conforme a especificação OpenSMILES) e quimicamente plausível, "
     "informando ainda em qual camada e posição a cadeia foi rejeitada.")

para("Principal função", bold=True)
para("Validar cadeias SMILES em camadas complementares: (1) análise léxica e sintática, "
     "por um parser LALR(1) próprio construído com a biblioteca SLY, cuja gramática foi "
     "adaptada do parser LL(1) de Richard Apodaca (depth-first.com) e da especificação "
     "OpenSMILES; (2) semântica de fechamento de anéis (OpenSMILES §3.4/§3.6), que exige "
     "estado ilimitado e por isso fica fora da gramática livre de contexto; e "
     "(3) validação química (valência e aromaticidade) sobre o grafo molecular, delegada "
     "a um motor químico maduro (pysmiles, por padrão) através de um adaptador neutro.")

para("Funções específicas relevantes oferecidas", bold=True)
bullet("análise léxica (com estados de bracket) e sintática de SMILES, com detecção "
       "precisa de erros e relato da camada e posição da falha.", "Parsing conforme OpenSMILES: ")
bullet("rejeita anéis não fechados, ligação de um átomo consigo mesmo, ligações "
       "duplicadas por fechamento de anel e ordens de ligação de anel incompatíveis.", "Semântica de anéis: ")
bullet("valida valência e aromaticidade delegando a pysmiles (padrão, permissivo) ou "
       "RDKit (estrito), sem nunca invocar o parser de strings do backend.", "Validação química: ")
bullet("validate_smiles devolve a tupla legada (bool, exceção); validate_smiles_detailed "
       "devolve um ValidationResult com camada (lex, grammar, ring_semantics, chemistry), "
       "mensagem e posição do caractere.", "Reporte por camadas: ")
bullet("uma tabela de divergências em relação à especificação OpenSMILES, congelada como "
       "teste de regressão (cada linha marcada com a seção da spec que a decide).", "Conformidade executável: ")

para("Usuários primordialmente contemplados", bold=True)
bullet("pesquisadores e engenheiros de IA generativa que precisam de um \"portão\" de "
       "validação para filtrar moléculas geradas por modelos generativos;")
bullet("químicos para farmacos e pesquisadores de quimioinformática que curam datasets de "
       "moléculas a partir de fontes heterogêneas;")
bullet("desenvolvedores que desejam reutilizar, total ou parcialmente, um parser LALR(1) "
       "de SMILES conforme OpenSMILES, desacoplado da validação química.")

para("Natureza do programa", bold=True)
para("Prova de conceito funcional e módulo em desenvolvimento, componente de uma "
     "iniciativa mais ampla de geração de moléculas com IA generativa. O parser sintático "
     "está maduro (gramática sem conflitos shift/reduce) e a validação química é delegada "
     "a um backend externo por um adaptador. O objetivo declarado é conformidade com a "
     "especificação OpenSMILES — o RDKit é usado como ponto de referência nas comparações, "
     "não como verdade absoluta.")

para("Ressalvas de uso", bold=True)
bullet("A validação química adota, por padrão, uma política deliberadamente permissiva: "
       "aceita espécies hipervalentes, radicais e íons livres que o RDKit rejeita; bloqueia "
       "apenas o quimicamente impossível (ex.: carbono pentavalente). Para a semântica "
       "conservadora do RDKit, selecione backend=\"rdkit\".")
bullet("Onde a especificação OpenSMILES diverge do RDKit, o validador segue a especificação. "
       "Algumas dessas divergências são intencionais e cobertas por testes (ver §3.1 e "
       "Apêndice A) — não devem ser \"corrigidas\" em direção ao RDKit.")
bullet("Aspectos estereoquímicos são reconhecidos (inclusive classes de quiralidade "
       "estendidas @TH/@AL/@SP/@TB/@OH), mas não influenciam a validação de valência e "
       "aromaticidade.")
bullet("Erros internos (backend ausente, entrada que não é string) são levantados como "
       "exceção — nunca reportados como \"molécula inválida\".")

# =========================================================================
# 2. VISÃO DE PROJETO (CENÁRIOS)
# =========================================================================
h("Visão de Projeto", 1)
para("Os cenários a seguir orientam o projeto, o uso e a evolução do programa, "
     "explicitando a intenção do criador e o comportamento esperado nas interações — "
     "inclusive nas que expõem limitações conhecidas.", italic=True)

h("Cenário Positivo 1 — Filtragem de moléculas geradas por IA", 2)
para("Marina desenvolve um modelo generativo que propõe novas moléculas candidatas a "
     "fármacos e recebe do modelo milhares de strings SMILES por lote. Antes de gastar "
     "tempo computacional em simulações caras, ela integra o Smiles Validator ao pipeline: "
     "para cada SMILES chama validate_smiles e mantém apenas as cadeias aceitas. As "
     "malformadas (parênteses abertos, anéis não fechados, átomos impossíveis) são "
     "rejeitadas com uma mensagem que indica a camada e a posição do erro, enquanto "
     "moléculas quimicamente reais — mesmo hipervalentes ou incomuns — passam. Marina "
     "obtém um lote limpo sem inspecionar manualmente milhares de strings.")

h("Cenário Positivo 2 — Curadoria de dataset para química medicinal", 2)
para("Otávio é químico medicinal e monta um dataset de compostos candidatos a partir de "
     "várias fontes públicas (ZINC, ChEMBL, extrações de literatura), reunindo centenas de "
     "milhares de strings SMILES de qualidade heterogênea. Ele passa cada SMILES por "
     "validate_smiles_detailed: as entradas malformadas são descartadas e ele registra a "
     "camada que as rejeitou (lex, grammar, ring_semantics ou chemistry), o que lhe permite "
     "auditar a origem dos problemas do dataset. Moléculas reais mas incomuns (hipervalentes, "
     "íons livres, radicais) são preservadas graças à política permissiva. Ao final, Otávio "
     "obtém um dataset limpo e confiável, sem ter descartado compostos legítimos que um "
     "validador conservador como o RDKit rejeitaria.")

h("Cenário Positivo 3 — Verificação de viabilidade de uma molécula projetada", 2)
para("A Dra. Sofia projeta, no papel, uma molécula que deveria exibir certas propriedades "
     "(por exemplo, um scaffold aromático com substituintes específicos para se ligar a um "
     "alvo). Antes de investir em síntese ou em simulações caras de propriedade, ela escreve "
     "as variantes candidatas em SMILES e as submete ao Smiles Validator para confirmar que "
     "representam moléculas estruturalmente possíveis — valências coerentes, anéis que fecham "
     "corretamente e aromaticidade válida. As variantes que violam a química básica são "
     "descartadas de imediato, e ela segue apenas com os esqueletos viáveis. Ressalva "
     "conhecida: o validador confirma a validade química e sintática da estrutura, não que a "
     "molécula terá as propriedades desejadas — isso exige ferramentas de predição a jusante; "
     "o validador apenas elimina cedo as ideias estruturalmente impossíveis.")

h("Cenário Negativo 1 — Divergência intencional com o RDKit (conformidade)", 2)
para("Bruno espera que o validador se comporte exatamente como o RDKit e submete "
     "C-1CCCCC=1, em que o símbolo de ligação do fechamento de anel é escrito nas duas "
     "pontas com ordens diferentes. O Smiles Validator rejeita a cadeia na camada de "
     "semântica de anéis, porque a OpenSMILES §3.4 exige que as duas pontas concordem; o "
     "RDKit, em contraste, silenciosamente deixa um dos símbolos vencer e aceita. Bruno "
     "também estranha que [C@@TH2](F)(Cl)(Br)I seja aceito aqui e recusado pelo RDKit. "
     "Ambos os casos são divergências intencionais, guiadas pela especificação e fixadas "
     "por testes dedicados: o alvo do projeto é a conformidade com o OpenSMILES, não a "
     "paridade com o RDKit. Se Bruno precisa da semântica do RDKit para a parte química, "
     "deve usar o backend \"rdkit\"; para a sintaxe, a referência é a especificação.")

h("Cenário Negativo 2 — Política de valência permissiva", 2)
para("Helena espera um validador químico estrito e submete N(C)(C)(C)C (nitrogênio neutro "
     "tetravalente) e alguns radicais. O Smiles Validator os aceita, porque a política "
     "padrão só rejeita a supervalência fisicamente impossível (ex.: carbono pentavalente) "
     "— a OpenSMILES não impõe um modelo de valência, e a saída de modelos generativos "
     "legitimamente contém radicais e espécies hipervalentes. Helena interpreta a aceitação "
     "como falha. Na verdade é uma decisão de projeto documentada: o portão foi concebido "
     "para não descartar moléculas reais. Ciente disso, Helena passa a usar "
     "backend=\"rdkit\" quando precisa da rigidez do sanitizador do RDKit, e mantém o padrão "
     "permissivo para triagem de moléculas geradas.")

# =========================================================================
# 3. DOCUMENTAÇÃO TÉCNICA
# =========================================================================
h("Documentação Técnica do Projeto", 1)
para("Esta seção destina-se a quem deseja reutilizar, total ou parcialmente, o programa. "
     "Inclui: (3.1) requisitos funcionais e não-funcionais; (3.2) arquitetura e modelo de "
     "dados; (3.3) modelo funcional; e (3.4) notas sobre o código.", italic=True)

# 3.1 Requisitos
h("3.1 Especificação de Requisitos", 2)
para("Requisitos funcionais (RF)", bold=True)
bullet("Expor um ponto de entrada único, validate_smiles(mol) -> (bool, Exceção|None), "
       "reexportado da raiz do pacote e de syntax.yacc (as duas rotas históricas não podem "
       "mais divergir).", "RF1 — ")
bullet("Oferecer validate_smiles_detailed(mol, backend) -> ValidationResult, com a camada "
       "que rejeitou (lex, grammar, ring_semantics, chemistry), mensagem, posição e "
       "tokens esperados.", "RF2 — ")
bullet("Aceitar a linguagem definida pela especificação OpenSMILES; rejeitar cadeias "
       "malformadas na camada correta.", "RF3 — ")
bullet("Aplicar a semântica de fechamento de anéis (§3.4/§3.6): anéis não fechados, "
       "self-ring-bond, ligação duplicada por anel e ordens de anel incompatíveis.", "RF4 — ")
bullet("Validar valência e aromaticidade via backend externo, sem nunca invocar o parser "
       "de strings do backend; permitir troca de política (pysmiles permissivo ↔ rdkit "
       "estrito) sem mudar a assinatura pública.", "RF5 — ")
bullet("Reportar verdicts de validação como valores; levantar exceção para erros internos "
       "(backend ausente, entrada não-string), que nunca devem ser confundidos com "
       "\"molécula inválida\".", "RF6 — ")

para("Requisitos não-funcionais (RNF)", bold=True)
bullet("Conformidade — as divergências conhecidas em relação à especificação OpenSMILES "
       "estão congeladas como testes de regressão (tests/conformance/test_opensmiles_gaps.py), "
       "cada uma marcada com a seção da spec que a decide.", "RNF1 — ")
bullet("Corretude estrutural — gramática LALR(1) sem conflitos shift/reduce (eram 26 antes "
       "da reescrita, agora 0).", "RNF2 — ")
bullet("Manutenibilidade — backend padrão em Python puro (pysmiles), passível de pin/patch; "
       "adaptador neutro fino.", "RNF3 — ")
bullet("Testabilidade — suíte com 118 testes passando (pytest), majoritariamente unitários "
       "escritos à mão (estados do lexer, semântica de anéis, adaptador, API de validação), "
       "mais a tabela de conformidade congelada (14) e um smoke test do dataset; inclui "
       "guardas de \"bypass\" que provam que o parser do backend nunca é chamado.", "RNF4 — ")
bullet("Portabilidade — instalável via uv/pip como pacote editável; datasets em Parquet "
       "(pyarrow).", "RNF5 — ")

para("Divergências documentadas em relação à BNF estrita / ao RDKit", bold=True)
bullet("bare atoms de um caractere fora do subconjunto orgânico; (E2) hcount multi-dígito "
       "([CH12]); (E3) faixas de dígitos de carga mais largas, com magnitude < 15.", "Extensões — (E1) ")
bullet("fechamento de anel não pode cruzar um ponto (C1.C1 é rejeitado); (R2) magnitude de "
       "carga ≥ 15 é rejeitada.", "Restrições — (R1) ")

# 3.2 Arquitetura
h("3.2 Arquitetura e Modelo de Dados", 2)
para("O sistema é um pipeline em camadas com fronteira nítida entre sintaxe/semântica de "
     "anéis (do projeto) e química (do backend), conforme a Figura 1. Cada camada "
     "corresponde a um tier reportável em caso de falha:")
doc.add_picture(_PIPELINE_PNG, width=Inches(4.4))
doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
cap = para("Figura 1 — Pipeline de validação do Smiles Validator, com os tiers de falha.", italic=True, size=9)
cap.alignment = WD_ALIGN_PARAGRAPH.CENTER

para("Módulos principais", bold=True)
bullet("o Parser LALR(1) — léxico com estados (SmilesLex para a cadeia principal, "
       "BracketLex para o interior de [...]) e gramática SLY, tratados como uma única etapa "
       "de análise sintática.", "src/syntax/ (lex.py, yacc.py, parser_manager.py) — ")
bullet("aplica a semântica de fechamento de anéis (§3.4/§3.6) e constrói o grafo canônico "
       "(âncora + pilha de ramificações), preservando a ordem de ligação também na cadeia "
       "principal.", "src/syntax/parser_manager.py — ")
bullet("modelo de dados químico — atomic.py (Atom/BracketAtom), structure.py "
       "(MolecularGraph, lista de adjacência), graph_builder.py (nós, arestas, anéis).", "src/chem/ — ")
bullet("tradução MolecularGraph → descritor neutro e drivers de backend (valid_pysmiles, "
       "valid_rdkit), bypassando o parser de strings do backend.", "src/chem/adapter.py — ")
bullet("montagem do pipeline e ponto de entrada público (validate_smiles, "
       "validate_smiles_detailed, ValidationResult, parse_smiles).", "src/validation.py — ")

para("Modelo de dados", bold=True)
bullet("ValidationResult (dataclass): valid (bool), tier, message, position, expected e "
       "error; __bool__ reflete valid e as_tuple() devolve a tupla legada.")
bullet("MolecularGraph: nós (símbolo, aromaticidade, carga, isótopo, hidrogênios "
       "explícitos, quiralidade) e arestas (ordem de ligação: '-', '=', '#', '$', ':' "
       "aromática — mapeadas para 1, 2, 3, 4 e 1,5 no descritor neutro).")

# 3.3 Modelo Funcional
h("3.3 Modelo Funcional", 2)
para("Fluxo de validação de uma cadeia SMILES:")
numbered("Análise sintática (Parser LALR(1)) — a string (com espaços nas bordas tolerados) "
         "é tokenizada por um lexer com estados: fora de colchetes vale o subconjunto "
         "orgânico e dígitos são fechamentos de anel; dentro de [...] vale a tabela "
         "periódica inteira, ':' é separador de classe atômica e sequências de dígitos "
         "viram um único token NUMBER (sem limite de largura). Erros de léxico/gramática "
         "geram falha nos tiers lex/grammar.")
numbered("Semântica de anéis — parser_manager.validate() verifica as regras §3.4/§3.6 "
         "(nenhum anel aberto etc.); violações geram RingSemanticsException (tier "
         "ring_semantics).")
numbered("Montagem do grafo — get_graph() devolve o MolecularGraph canônico já construído "
         "durante as ações semânticas.")
numbered("Tradução — to_descriptor converte o grafo no descritor neutro (atoms, bonds).")
numbered("Química — o driver (pysmiles por padrão) executa percepção de anéis, "
         "kekulização/aromaticidade e o modelo de valência multi-valente; falhas do backend "
         "caem no tier chemistry.")
numbered("Retorno — ValidationResult (truthy sse válido); validate_smiles devolve a tupla "
         "legada via as_tuple(). Erros internos são levantados, não reportados como inválido.")

para("Regras químicas de referência", bold=True)
bullet("Regra de Hückel (aromaticidade): um ciclo planar conjugado é aromático se possui "
       "4n+2 elétrons π (2, 6, 10, 14, ...). Percepção de anéis e kekulização são do backend.")
bullet("Política de valência (permissiva): aceita radicais e espécies de octeto expandido / "
       "hipervalentes; rejeita apenas supervalência impossível. Opt-in para rigidez do "
       "RDKit via backend=\"rdkit\".")

# 3.4 Sobre o código
h("3.4 Sobre o Código", 2)
bullet("Python 3.10+; parser gerado com o toolkit SLY (LALR(1)); grafos com NetworkX; "
       "dados em Parquet (pandas/PyArrow); química com pysmiles e, opcionalmente, RDKit.", "Linguagem e bibliotecas: ")
bullet("lexer com estados (SmilesLex/BracketLex) e troça de estado ao encontrar '['/']'; "
       "tokens NUMBER (sem limite) e CHIRAL (@, @@ e classes @TH/@AL/@SP/@TB/@OH); "
       "combinações de campos opcionais do bracket geradas por generate_combinations; "
       "tabela periódica carregada de periodic-table-lookup.json.", "Técnica: ")
bullet("docstrings em estilo Google renderizadas via mkdocstrings/MkDocs Material; "
       "comentários em linha reservados a decisões não óbvias (ex.: por que a semântica de "
       "anéis fica fora da gramática, ou por que a valência é permissiva).", "Estratégia de comentários: ")
bullet("conformidade executável — tests/conformance/test_opensmiles_gaps.py congela a tabela "
       "de divergências em relação à especificação OpenSMILES como teste de regressão, cada "
       "linha marcada com a seção da spec que a decide.", "Conformidade: ")
bullet("adaptador com contrato estável; testes de \"bypass\" (monkeypatch de "
       "read_smiles/MolFromSmiles) garantem que o parser de strings do backend jamais "
       "entra no caminho de execução.", "Extensibilidade: ")
bullet("pytest (118 testes); documentação em https://luna-v0.github.io/smiles/.", "Testes e docs: ")

# =========================================================================
# 4. MANUAL DE UTILIZAÇÃO
# =========================================================================
h("Manual de Utilização para Usuários Contemplados", 1)

h("Instalação", 2)
para("Requer o gerenciador uv (astral-sh/uv). No diretório do projeto:")
code_block([
    "uv venv .venv",
    "source .venv/bin/activate",
    "uv pip install -e .",
])

h("Tarefa A — Validar uma cadeia SMILES (uso mais comum)", 2)
para("Guia de instruções:", bold=True)
numbered("Importe a função da raiz do pacote: from src import validate_smiles "
         "(também disponível, por compatibilidade, em syntax.yacc).")
numbered("Chame validate_smiles(\"<sua_cadeia_SMILES>\").")
numbered("Inspecione o retorno: uma tupla (bool, Exceção|None).")
code_block([
    ">>> from src import validate_smiles",
    ">>> validate_smiles(\"c1ccccc1\")      # benzeno",
    "(True, None)",
    ">>> validate_smiles(\"C1CCC\")         # anel não fechado",
    "(False, RingSemanticsException(...))",
])
para("Exceções ou potenciais problemas:", bold=True)
bullet("Espaços/quebras de linha nas bordas são tolerados (a entrada é normalizada com "
       "strip); caracteres ilegais no meio da cadeia caem no tier lex.", "C1: ")
bullet("Se a entrada não for uma string, ou o backend químico não estiver instalado, a "
       "função levanta exceção (TypeError/ImportError) — isso não é significa que a "
       "molecula é \"inválida\".", "C2: ")

h("Tarefa B — Diagnóstico detalhado (camada e posição da falha)", 2)
para("Guia de instruções:", bold=True)
numbered("Importe from src import validate_smiles_detailed.")
numbered("Chame-a e leia os campos do ValidationResult (valid, tier, message, position).")
code_block([
    ">>> from src import validate_smiles_detailed",
    ">>> r = validate_smiles_detailed(\"C1CCC\")",
    ">>> r.valid, r.tier, r.message",
    "(False, 'ring_semantics', 'Unclosed ring numbers: 1')",
    ">>> r.as_tuple()   # contrato legado (bool, exceção)",
])
para("Use esta tarefa quando precisar auditar por que uma cadeia foi rejeitada (ex.: "
     "curadoria de dataset), distinguindo erro de léxico, gramática, anel ou química.", italic=True)

h("Tarefa C — Escolher a política de validação química", 2)
para("Guia de instruções:", bold=True)
numbered("Passe backend=\"pysmiles\" (permissivo, padrão) ou backend=\"rdkit\" (estrito) a "
         "validate_smiles_detailed, ou instancie ChemistryValidator(backend=...).")
para("Em contextos diferentes: para triagem de moléculas geradas por IA, prefira "
     "\"pysmiles\" (aceita moléculas reais incomuns); para máxima concordância com o RDKit "
     "como referência, prefira \"rdkit\".", italic=True)

h("Tarefa D — Reutilizar apenas o parser (sem química)", 2)
para("Guia de instruções:", bold=True)
numbered("Importe from src import parse_smiles.")
numbered("Chame parse_smiles(\"<cadeia>\"); recebe (is_valid, graph, exceção) com a sintaxe "
         "e a semântica de anéis validadas, mas sem validação química.")
para("Exceções ou potenciais problemas:", bold=True)
bullet("Casos aromáticos fundidos muito complexos continuam sendo delegados ao backend na "
       "etapa química; parse_smiles apenas garante sintaxe + anéis.", "C1: ")

h("Tarefa E — Executar os testes", 2)
para("Guia de instruções:", bold=True)
numbered("Com o ambiente virtual ativo, execute pytest na raiz do projeto.")
code_block(["pytest    # 118 testes"])

# =========================================================================
# Referências
# =========================================================================
h("Referências", 1)
bullet("OpenSMILES Specification — https://opensmiles.org/opensmiles.html")
bullet("R. Apodaca, \"SMILES Formal Grammar\" — https://depth-first.com/articles/2020/04/20/smiles-formal-grammar/")
bullet("SLY (David Beazley) — https://github.com/dabeaz/sly")
bullet("Bowserinator/Periodic-Table-JSON — dados da tabela periódica")
bullet("Documentação do projeto — https://luna-v0.github.io/smiles/")
bullet("J. Carroll (1999), \"Five Reasons for Scenario-Based Design\" — base metodológica dos cenários")

# =========================================================================
# Apêndice A — Gramática
# =========================================================================
doc.add_page_break()
h("Apêndice A — Gramática do Parser LALR(1)", 1)
para("A gramática livre de contexto reconhecida pelo parser LALR(1) (SLY) é reproduzida "
     "abaixo, transcrita da saída do gerador (parser.out). Notação: \"→\" produção, \"|\" "
     "alternativas, sufixo \"?\" opcional, aspas simples para terminais literais. Os "
     "terminais NUMBER (sequência de dígitos, sem limite) e CHIRAL (designador de "
     "quiralidade) são emitidos pelo lexer no estado de bracket.", italic=True)

code_block([
    "line             → atom",
    "                 | atom chain_branch",
    "",
    "chain_branch     → chains",
    "                 | branch",
    "                 | branch chain_branch",
    "                 | chains chain_branch",
    "",
    "chains           → chain",
    "                 | chain chains",
    "",
    "chain            → rnum",
    "                 | atom",
    "                 | bond '%' digit digit",
    "                 | bond digit",
    "                 | bond atom",
    "                 | dot_proxy",
    "",
    "dot_proxy        → '.' atom",
    "",
    "atom             → bracket_atom",
    "                 | symbol",
    "",
    "bracket_atom     → '[' internal_bracket ']'",
    "internal_bracket → isotope? symbol chiral? hcount? charge? mol_map?",
    "                   (symbol obrigatório; demais opcionais — expandido em",
    "                    32 regras explícitas, 10–41, no parser.out)",
    "",
    "symbol           → '*' | 'H' | semi_symbol",
    "",
    "branch           → open_paren inner_branch close_paren",
    "open_paren       → '('",
    "close_paren      → ')'",
    "",
    "inner_branch     → line inner_branch",
    "                 | branch_bond line inner_branch",
    "                 | line",
    "                 | branch_bond line",
    "",
    "branch_bond      → bond_dot",
    "bond_dot         → '.' | bond",
    "bond             → ':' | '-' | semi_bond_rule",
    "semi_bond_rule   → semi_bond",
    "",
    "rnum             → '%' digit digit | digit",
    "isotope          → NUMBER",
    "hcount           → 'H' NUMBER | 'H'",
    "charge           → '+' | '-' | '+' '+' | '-' '-'",
    "                 | '+' fifteen | '-' fifteen",
    "mol_map          → ':' NUMBER",
    "chiral           → CHIRAL",
    "fifteen          → NUMBER",
])

para("Terminais", bold=True)
bullet("símbolos de elemento da tabela periódica e símbolos aromáticos "
       "(b, c, n, o, p, s, se, as).", "semi_symbol — ")
bullet("ligações não-simples fora do bracket: '=' (dupla), '#' (tripla), '$' (quádrupla), "
       "'/' e '\\' (estéreo). A ligação aromática ':' também é aceita como bond.", "semi_bond — ")
bullet("sequência de um ou mais dígitos, sem limite de largura (isótopo, hcount, carga, "
       "classe atômica), emitida dentro de [...].", "NUMBER — ")
bullet("designador de quiralidade: '@', '@@' ou classe @TH/@AL/@SP/@TB/@OH com índice "
       "(faixa validada semanticamente).", "CHIRAL — ")
bullet("dígito 0–9 usado em números de anel fora do bracket.", "digit — ")
bullet("'.' '@' '-' '+' ':' '%' 'H' '(' ')' '[' ']' '*'.", "Literais — ")

para("Observações", bold=True)
bullet("O símbolo inicial é line (o parser.out o expõe como S' → line).")
bullet("A gramática não possui conflitos shift/reduce (0), ante os 26 da versão anterior; "
       "a semântica de fechamento de anéis, que exige estado ilimitado, fica fora da "
       "gramática, no ParserManager (OpenSMILES §3.4/§3.6).")
bullet("As 32 produções de internal_bracket são geradas automaticamente por "
       "generate_combinations a partir do padrão isotope? symbol chiral? hcount? charge? "
       "mol_map?, sempre com symbol presente.")

doc.save(OUT)
print("OK:", OUT)
