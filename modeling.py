import marimo

__generated_with = "0.23.4"
app = marimo.App()


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _():
    import os
    print(os.system("python -V"))
    return (os,)


@app.cell
def _(os):
    print(os.environ.get('ILOG_LICENSE_FILE'))
    print(os.environ.get('CPLEX_STUDIO_DIR'))
    # Check environment variables
    for key, value in os.environ.items():
        if 'cplex' in key.lower() or 'ilog' in key.lower() or 'ibm' in key.lower():
    # Print all CPLEX/ILOG related environment variables
            print(f'{key}={value}')
    return


@app.cell
def _(os):
    import cplex
    c = cplex.Cplex()
    print(c.get_version())
    # Some installations store info here
    print(os.path.dirname(cplex.__file__))
    return


@app.cell
def _():
    import subprocess
    _result = subprocess.run(['find', '/opt/env/modelseed_cplex', '-name', 'access.ilm'], capture_output=True, text=True)
    # Search within the virtual environment
    print(_result.stdout)
    _result = subprocess.run(['find', '/opt/env/modelseed_cplex', '-name', '*.ilm'], capture_output=True, text=True)
    # Also check for any .ilm files
    print(_result.stdout)
    return (subprocess,)


@app.cell
def _(subprocess):
    _result = subprocess.run(['find', '/opt', '-name', 'access.ilm'], capture_output=True, text=True, timeout=120)
    print(_result.stdout)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # create models from ASVset genome objects
    """)
    return


@app.cell
def _(missing_IDs, missing_iterativeIDs, os):
    from load_genomes_to_modelseedpy import GenomeLoader
    from cobra.io import write_sbml_model
    from tqdm import tqdm
    loader = GenomeLoader()

    def write_model(gID):
        model = loader.build_model(gID, printing=False)
        write_sbml_model(model, f'models/{gID}.xml')
    _args = [gID for gID in loader.list_available_genomes() if not os.path.exists(f'models/{gID}.xml')]
    print(len(_args), 'will be processed')
    _parallelize = False
    if _parallelize:
        from multiprocess import Pool
        from os import cpu_count
        _cpus = int(cpu_count() / 8)
        print(f'{_cpus} cores are being used.  The first argument is {_args[0]}')
        _pool = Pool(_cpus)
        _outputs = _pool.map(write_model, _args)
        print(f'missing iterativeIDs', missing_iterativeIDs)
        print(f'missing IDs', missing_IDs)
    else:
        for _arg in _args:
            write_model(_arg)
    return Pool, cpu_count, tqdm, write_sbml_model


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Gapfill the models in minimal media
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### ATP Correction
    """)
    return


@app.cell
def _():
    from cobra.io import read_sbml_model
    from glob import glob

    member_models = {}
    num_models = 5
    for model in glob('models/*.xml')[:num_models]:
        ID = model.split("/")[0].split(".")[0]
        member_models[ID] = read_sbml_model(model)
    return glob, member_models, read_sbml_model


@app.cell
def _(member_models):
    from modelseedpy.core.mstemplate import MSTemplateBuilder
    from modelseedpy import MSATPCorrection
    from json import load
    with open('../ModelSEEDpy/modelseedpy/data/templates/template_core.json') as _fh:
        template_core = MSTemplateBuilder.from_dict(load(_fh)).build()
    ATPcorrection_tests = {}
    for ID_1, model_1 in member_models.items():
        atp_correction = MSATPCorrection(model_1, template_core, compartment='c0', atp_hydrolysis_id='ATPM_c0', load_default_medias=True)
        print('correction done')
        media_eval = atp_correction.evaluate_growth_media()
        print('evaluation done')
        atp_correction.determine_growth_media()
        print('growth media done')
        atp_correction.apply_growth_media_gapfilling()
        print('applying growth media done')
        atp_correction.expand_model_to_genome_scale()
        print('expanding model done')
        ATPcorrection_tests[ID_1] = atp_correction.build_tests()
    return ATPcorrection_tests, MSATPCorrection, MSTemplateBuilder, load


@app.cell
def _(load, util):
    # magic command not supported in marimo; please file an issue to add support
    # %run util.py
    from numpy import mean
    microbiome_path = '../MicrobiomeNotebooks'
    metabolites = load(open(f'{microbiome_path}/digestor/metabolites.json', 'r'))
    # microbiome_path = "/Users/andrewfreiburger/Documents/Research/MicrobiomeNotebooks/digestor"
    metaboliteIDs = list(metabolites.keys())
    media_nar = 207617
    metabolites = util.msrecon.get_media(f'{media_nar}/Wolfe').mediacompounds
    auxo_media = util.msrecon.get_media(f'{media_nar}/AuxoMedia')
    anaerobicPyruvate = util.msrecon.get_media(f'{media_nar}/PyruateMinimalAnaerobic')
    metaboliteIDs = [cpd['id'] for cpd in metabolites] + ['cpd00011', 'cpd11640', 'cpd01024']
    uptake_phenoset = util.create_phenotypeset_from_compounds(metaboliteIDs, base_media=auxo_media, base_uptake=0, base_excretion=1000, global_atom_limits={}, type='uptake')
    excretion_phenoset = util.create_phenotypeset_from_compounds(metaboliteIDs, base_media=auxo_media, base_uptake=0, base_excretion=1000, global_atom_limits={}, type='excretion')
    # metabolites
    agrowth_phenoset = util.create_phenotypeset_from_compounds(metaboliteIDs, base_media=anaerobicPyruvate, base_uptake=0, base_excretion=1000, global_atom_limits={}, type='agrowth')  # experimental media + CO2,H2,CH4
    #TODO:  Why is growth only performed for the non-base media and uptake/excretion only performed for the base media?
    # growth_phenoset = util.create_phenotypeset_from_compounds(
    #     metabolites, base_media=gmm_base_media, base_uptake=0, base_excretion=1000, global_atom_limits={}, type="growth")
    phenotypes = {'uptake': uptake_phenoset, 'excretion': excretion_phenoset, 'agrowth': agrowth_phenoset}
    return anaerobicPyruvate, auxo_media


@app.cell
def _(
    ATPcorrection_tests,
    MSTemplateBuilder,
    anaerobicPyruvate,
    auxo_media,
    gID,
    genome,
    json,
    member_models,
    write_sbml_model,
):
    from modelseedpy.core.msmodel import get_reaction_constraints_from_direction
    from modelseedpy import MSBuilder
    from modelseedpy import MSGapfill
    from modelseedpy.core.mspredict import MSPredict

    def _integrate_solution(template, model, gap_fill_solution):
        added_reactions = []
        for rxn_id, (lb, ub) in gap_fill_solution.items():
            template_reaction = template.reactions.get_by_id(rxn_id)
            model_reaction = template_reaction.to_reaction(model)
            model_reaction.lower_bound = lb
            model_reaction.upper_bound = ub
            _str = model_reaction.build_reaction_string(True)
            added_reactions.append(model_reaction)
        model.add_reactions(added_reactions)
        add_exchanges = MSBuilder.add_exchanges_to_model(model)
        return (added_reactions, add_exchanges)
    with open('~/repos/ModelSEEDpy/modelseedpy/data/templates/template_gram_neg.json') as _fh:
        template_gramneg = MSTemplateBuilder.from_dict(json.load(_fh)).build()
    for ID_2, model_2 in member_models:
        predict = MSPredict()
        genome_klass = predict.predict(genome)
        gapfill = MSGapfill(model_2, default_gapfill_templates=[template_gramneg], test_conditions=ATPcorrection_tests[ID_2], default_target='bio1')
        for media in [auxo_media, anaerobicPyruvate]:
            gapfill_res = gapfill.run_gapfilling(auxo_media)
            gap_sol = {}
            for rxn_id, d in gapfill_res['new'].items():
                if rxn_id[:-1] in template_gramneg.reactions:
                    gap_sol[rxn_id[:-1]] = get_reaction_constraints_from_direction(d)
            print(gap_sol)
            _integrate_solution(template_gramneg, model_2, gap_sol)
        write_sbml_model(model_2, f'models/{gID}.xml')
    return (MSPredict,)


@app.cell
def _(MSPredict, glob, os, tqdm):
    # Cell 1: Imports and setup                                                                                                                
    import json
    from collections import Counter
    from modelseedpy.core.msgenome import MSGenome, MSFeature
    from gapfill_models import get_template_for_genome_class, create_auxo_media, create_pyruvate_minimal_anaerobic, gapfill_model, list_models
    models_dir = 'models'
    genome_objects_dir = 'genome_objects'
    output_dir = 'gapfilled_models2'
    predictor = MSPredict()

    def load_genome_from_genome_object(genome_id):
        """Load MSGenome with RAST annotations from genome_objects JSON."""
        obj_path = os.path.join(genome_objects_dir, f'{genome_id}.json')
        if not os.path.exists(obj_path):
            return None
        genome = MSGenome()
        with open(obj_path) as f:
            data = json.load(f)
        genome.id = data.get('id', genome_id)
        features = data.get('features', [])
        for feat in features:
    # Initialize MSPredict classifier                                                                                                          
            if not isinstance(feat, dict):
                continue
            if feat.get('type') == 'CDS':
                feat_id = feat.get('id', '')
                ms_feat = MSFeature(feat_id, '')
                ontology = feat.get('ontology_terms', {})
                if 'RAST' in ontology and ontology['RAST']:
                    for term in ontology['RAST']:
                        ms_feat.add_ontology_term('RAST', term)
                    genome.features.append(ms_feat)
                elif feat.get('functions'):
                    for func in feat['functions']:
                        ms_feat.add_ontology_term('RAST', func)
                    genome.features.append(ms_feat)
        return genome
    genome_ids = [os.path.splitext(os.path.basename(f))[0] for f in glob(f'{genome_objects_dir}/*.json')]
    genome_classifications = {}  # Check for CDS type and has functions or ontology_terms                                                                           
    print('classifying the genomes')
    for gid in tqdm(genome_ids):
        genome = load_genome_from_genome_object(gid)
        if genome and len(genome.features) > 0:  # Use ontology_terms if available, otherwise use functions                                                                     
            try:
                genome_classifications[gid] = predictor.predict(genome).value
            except Exception as e:
                print(f'Could not classify {gid}: {e}')
    json.dump(genome_classifications, open('member_classifications.json', 'w'), indent=3)
    counts = Counter(list(genome_classifications.values())).most_common()
    # Get genome IDs from genome_objects directory                                                                                             
    print(f'Classified {len(genome_classifications)} genomes: {counts}')
    return (
        create_auxo_media,
        create_pyruvate_minimal_anaerobic,
        gapfill_model,
        genome,
        genome_classifications,
        json,
        list_models,
        models_dir,
        output_dir,
        predictor,
    )


@app.cell
def _(genome_classifications):
    list(genome_classifications.items())[:2]
    return


@app.cell
def _(json, list_models, models_dir, os, tqdm):
    genome_classifications_1 = json.load(open('member_classifications.json', 'r'))
    print('ID mapping')
    model_genome_mapping = {}
    for model_path in tqdm(list_models(models_dir)):
        _model_id = os.path.splitext(os.path.basename(model_path))[0]
        if _model_id in genome_classifications_1:
            model_genome_mapping[_model_id] = _model_id
        else:
            base_id = '.'.join(_model_id.split('.')[:2])
            if base_id in genome_classifications_1:
                model_genome_mapping[_model_id] = base_id
    print(f'Mapped {len(model_genome_mapping)} models to genomes')
    return genome_classifications_1, model_genome_mapping


@app.cell
def _(
    MSATPCorrection,
    Pool,
    cpu_count,
    create_auxo_media,
    create_pyruvate_minimal_anaerobic,
    gapfill_model,
    genome_classifications_1,
    list_models,
    model_genome_mapping,
    models_dir,
    os,
    output_dir,
    predictor,
    read_sbml_model,
    write_sbml_model,
):
    print('Gapfilling')
    from modelseedpy.core.mspredict import MSGenomeClass
    media_list = [create_auxo_media(), create_pyruvate_minimal_anaerobic()]
    atp_correction_tests = {}
    results = {}
    genome_class_map = {'Gram Positive': MSGenomeClass.P, 'Gram Negative': MSGenomeClass.N, 'Cyano': MSGenomeClass.C, 'Archaea': MSGenomeClass.A}

    def gapfilling(model_path):
        _model_id = os.path.splitext(os.path.basename(model_path))[0]
        genome_id = model_genome_mapping.get(_model_id)
        genome_class_str = genome_classifications_1.get(genome_id) if genome_id else None
        genome_class = genome_class_map.get(genome_class_str) if genome_class_str else None
        core_template = None
        genome_template = None
        if genome_class:
            core_template, genome_template = predictor.auto_select_template(genome_class)
        output_path = gapfill_model(model_path=model_path, media_list=media_list, output_dir=output_dir, gapfilling_mode='Sequential', atp_safe=True, genome_class=genome_class)
        _result = None
        if output_path:
            model = read_sbml_model(output_path)
            if core_template:
                atp_correction = MSATPCorrection(model, core_template, compartment='c0', atp_hydrolysis_id='ATPM_c0', load_default_medias=True)
                atp_correction.evaluate_growth_media()
                atp_correction.determine_growth_media()
                atp_correction.apply_growth_media_gapfilling()
                atp_correction.expand_model_to_genome_scale()
                tests = atp_correction.build_tests()
                for t in tests:
                    print(t['media'].id, t['threshold'], atp_correction.media_gapfill_stats[t['media']])
                write_sbml_model(model, output_path)
                _result = (_model_id, output_path, tests)
            else:
                _result = (_model_id, output_path, None)
            print(f'Completed: {_model_id}')
        else:
            _result = (_model_id, None, None)
        return _result
    _args = list_models(models_dir)
    print(len(_args), 'will be processed')
    _parallelize = True
    if _parallelize:
        _cpus = min(len(_args), int(cpu_count() / 8))
        print(f'{_cpus} cores are being used.  The first argument is {_args[0]}')
        _pool = Pool(_cpus)
        _outputs = _pool.map(gapfilling, _args)
        _pool.close()
        _pool.join()
        for _model_id, output_path, tests in _outputs:
            results[_model_id] = output_path
            if tests:
                atp_correction_tests[_model_id] = tests
    else:
        for _arg in _args:
            _model_id, output_path, tests = gapfilling(_arg)
            results[_model_id] = output_path
            if tests:
                atp_correction_tests[_model_id] = tests
    print(f'Completed: {sum((1 for v in results.values() if v))}/{len(results)} successful')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Define the metabolite-centered phenotypes
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # create community models from the member models and abundances
    """)
    return


@app.cell
def _(
    Pool,
    cpu_count,
    load,
    missing_IDs,
    missing_iterativeIDs,
    os,
    read_sbml_model,
    write_sbml_model,
):
    from pandas import read_csv, Series
    from mscommunity import build_from_species_models
    iterativeIDs = load(open('model_inputs/iterativeIDs.json', 'r'))
    ASV_sets = load(open('modeling_files/ASV_sets.json', 'r'))
    redunant_asvs = {v: k for k, vs in ASV_sets.items() for v in vs}
    _abundances = load(open('modeling_files/ASVset_abundances.json', 'r'))
    _abundances = {sample: {iterativeIDs.get(ASV): abund for ASV, abund in content.items()} for sample, content in _abundances.items()}
    abundances_new = {sample: {iterativeIDs.get(ASV): {'abundance': abund} for ASV, abund in content.items()} for sample, content in _abundances.items()}

    def create_comm(item):
        sample, _abundances = item
        models = []
        for ASVset, abun in _abundances.items():
            model_path = f'gapfilled_models/{ASVset}_gf.xml'
            if not os.path.exists(model_path):
                continue
            model = read_sbml_model(model_path)
            model.id = ASVset
            models.append(model)
        modelID = f'{sample}_comm'
        modelName = f'{sample}_comm'
        print(len(models))
        comm_model = build_from_species_models(models, modelID, modelName, printing=True)
        write_sbml_model(comm_model, f'comm_models/{modelID}.xml')
    _args = [item for item in _abundances.items() if not os.path.exists(f'comm_models/{item[0]}_comm.xml')]
    print(len(_args), 'will be processed')
    _parallelize = False
    if _parallelize:
        _cpus = min(len(_args), int(cpu_count() / 8))
        print(f'{_cpus} cores are being used.  The first argument is {_args[0]}')
        _pool = Pool(_cpus)
        _outputs = _pool.map(create_comm, _args)
        print(f'missing iterativeIDs', missing_iterativeIDs)
        print(f'missing IDs', missing_IDs)
    else:
        for _arg in _args:
            create_comm(_arg)
    return (read_csv,)


@app.cell
def _(load):
    _abundances = load(open('modeling_files/ASVset_abundances.json', 'r'))
    list(list(_abundances.values())[0].values())[0]
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Test that the models are operational
    """)
    return


@app.cell
def _(read_sbml_model):
    # for path in glob("models/F34_comm.xml"):
    # for path in glob("comm_models/G12_comm.xml"):
    # break
    model_3 = read_sbml_model('comm_models/P12_comm.xml')  # print(path.split("/")[-1])
    return (model_3,)


@app.cell
def _(display, model_3):
    print(model_3.objective)
    display(model_3.summary())
    return


@app.cell
def _(model_3):
    model_3.reactions.bio1.reaction
    return


@app.cell
def _(model_3):
    bioReactions = []
    for rxn in model_3.reactions:
        if 'bio' in rxn.id:
            bioReactions.append(rxn.id)
    bioReactions = sorted(bioReactions)
    print(len(bioReactions), bioReactions)
    return


@app.cell
def _(model_3):
    bioCompounds = []
    for cpd in model_3.metabolites:
        if 'cpd11416' in cpd.id:
            bioCompounds.append(cpd.id)
    bioCompounds = sorted(bioCompounds)
    print(len(bioCompounds), bioCompounds)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Community modeling
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Load the community models
    """)
    return


@app.cell
def _(glob, read_sbml_model, tqdm):
    # magic command not supported in marimo; please file an issue to add support
    # %%time
    def load_models():
        models = {}
        num_models = 1
        modelID = 'P12_comm'
        for model_path in tqdm(glob(f'comm_models/{modelID}.xml')[:num_models]):
            modelID = model_path.split('/')[-1].replace('.xml', '')
            print('loading', modelID)
            model = read_sbml_model(model_path)
            model.id = modelID
            models[modelID] = model
        return models
    models = load_models()
    return load_models, models


@app.cell
def _(modelID, models):
    type(models[modelID].solver)
    return


@app.cell
def _(display, modelID, models):
    # magic command not supported in marimo; please file an issue to add support
    # %%time
    display(models)
    print(models[modelID].objective)
    display(models[modelID].summary())
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Constraints
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Process data
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    #### convert measurements into MSIDs
    """)
    return


@app.cell
def _(read_csv):
    from pandas import set_option
    from json import dump
    set_option('display.max_columns', None)
    mapping = {'Media Acetate (mol/L)': 'cpd00029_in', 'Media Propionate (mol/L)': 'cpd00141_in', 'Media Butyrate (mol/L)': 'cpd00211_in', 'Waste Effluent Acetate (mol/L)': 'cpd00029_out', 'Waste Effluent Propionate (mol/L)': 'cpd00141_out', 'Waste Effluent Butyrate (mol/L)': 'cpd00211_out', 'Approximate Total Influent C (mol/min) based on Media recipe': 'media_carbon_in', 'Measured Total Influent C (mol/min)': 'carbon_in', 'Gas Composition (% CH4)': '%cpd01024', 'Gas Composition (% H2)': '%cpd11640', 'Gas Composition (% CO2)': '%cpd00011', 'CH4 breakthrough (mol/min)': 'cpd01024_out', 'H2 breakthrough (mol/min)': 'cpd11640_out', 'CO2 breakthrough (mol/min)': 'cpd00011_out', 'Accounted C (%) based on MT sensors': 'accounted_C'}
    data = read_csv('model_inputs/measurements/Summary_interpolated.csv')
    data['Media Total COD (mg/L) Hach'] = data['Media Total COD (mg/L) Hach'] / 32000
    data['carbon_out'] = data['Waste Effluent Dissolved Inorganic Carbon (mol C/L)'] + data['Waste Effluent Dissolved Organic Carbon (mol C/L)'] + data['Effluent C as measured (mol/min)']
    data.rename(columns=mapping, inplace=True)
    data = data[data['Biomass Sample ID'].notna()].set_index('Biomass Sample ID')
    data_df = data.T.to_dict()
    dump(data_df, open('model_inputs/measurements/Summer_interpolated_MSDB.json', 'w'), indent=2)  # "Media Total COD (mg/L) Hach": "media_cpd00007",
    # print(Counter(data.columns))
    data
    return (data_df,)


@app.cell
def _(data_df):
    data_df
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### loading the media
    """)
    return


@app.cell
def _(load):
    # import os
    # # os.environ["HOME"] = "~/repos/cobrakbase"
    # import cobrakbase
    # with open("/home/afreiburger/.kbase/token", 'r') as token:
    #     kbase_obj = cobrakbase.KBaseAPI(token.readline())
    # media_ws = 207617
    # media = kbase_obj.get_from_ws("Wolfe",207617)
    # from modelseedpy.core.fbahelper import FBAHelper
    # from json import dump
    # media_js = FBAHelper.convert_kbase_media(media)
    # dump(media_js, open("model_inputs/Wolfe.json", 'w'), indent=2)
    media_1 = load(open('model_inputs/Wolfe.json', 'r'))
    return (media_1,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ###
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Adding the media and element constraints
    """)
    return


@app.cell
def _(display, load_models, media_1, tqdm):
    # magic command not supported in marimo; please file an issue to add support
    # %%time
    from mscommunity import MSCommunity
    models_1 = load_models()
    for ID_3, model_4 in tqdm(models_1.items()):
        sample = ID_3.split('_')[0]
        models_1[ID_3] = MSCommunity(model_4)
        display(models_1[ID_3].util.model.summary())
        models_1[ID_3].util.add_medium = media_1  #, climit=data_df[sample]["carbon_in"]*1000, o2limit=0)
    return ID_3, models_1


@app.cell
def _(ID_3, display, models_1):
    # magic command not supported in marimo; please file an issue to add support
    # %%time
    # Without Carbon and Oxygen constraints
    display(models_1[ID_3].util.model.summary())
    return


@app.cell
def _(ID_3, display, models_1):
    # magic command not supported in marimo; please file an issue to add support
    # %%time
    # With Carbon and Oxygen constraints
    display(models_1[ID_3].util.model.summary())
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Objective
    """)
    return


@app.cell
def _(models_1):
    # magic command not supported in marimo; please file an issue to add support
    # %%time
    for ID_4, model_5 in models_1.items():
        print(ID_4)
        print(model_5.util.model.objective.expression)
        print(model_5.run_fba())
    return


if __name__ == "__main__":
    app.run()
