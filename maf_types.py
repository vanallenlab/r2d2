class MafTypes(object):
    dna_normal = 'dna_normal'
    dna_tumor = 'dna_tumor'
    rna_normal = 'rna_normal'
    rna_tumor = 'rna_tumor'

    @classmethod
    def all(cls):
        return [MafTypes.dna_normal,
                MafTypes.dna_tumor,
                MafTypes.rna_normal,
                MafTypes.rna_tumor]


class AnalysisTypes(object):
    all_inputs = 'all_inputs'
    no_rna_normal = 'no_rna_normal'
    dna_only = 'dna_only'
    normal_only = 'normal_only'
    tumor_only = 'tumor_only'


ANALYSIS_SAMPLES = {
    AnalysisTypes.all_inputs: [MafTypes.dna_normal, MafTypes.dna_tumor, MafTypes.rna_normal, MafTypes.rna_tumor],
    AnalysisTypes.no_rna_normal: [MafTypes.dna_normal, MafTypes.dna_tumor, MafTypes.rna_tumor],
    AnalysisTypes.dna_only: [MafTypes.dna_normal, MafTypes.dna_tumor],
    AnalysisTypes.normal_only: [MafTypes.dna_normal, MafTypes.rna_normal],
    AnalysisTypes.tumor_only: [MafTypes.dna_tumor, MafTypes.rna_tumor]
}