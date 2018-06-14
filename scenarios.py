import ConfigParser

SCENARIOS_CONFIG = 'scenarios.ini'

class Condition(object):
    conditions_types = ['=', '<', '<=', '>', '>=', '<>']

    def __init__(self, condition_str):
        self._clauses = []

        for clause in condition_str.split('|'):
            tokens = clause.strip().split(' ')
            tokens = [tok.strip() for tok in tokens]

            if len(tokens) == 2:
                threshold = float(tokens[1])
            elif len(tokens) == 3:
                threshold = (float(tokens[1]), float(tokens[2]))
            else:
                raise ValueError('Malformed conditions string - incorrect length (%s)' % condition_str)

            if tokens[0] not in Condition.conditions_types:
                raise ValueError('Malformed conditions string - incorrect type (%s)' % condition_str)

            self._clauses.append({'type': tokens[0], 'threshold': threshold})

    def test(self, value):
        value = float(value)
        for clause in self._clauses:
            if (
                (clause['type'] == '<' and value < clause['threshold']) or
                (clause['type'] == '<=' and value <= clause['threshold']) or
                (clause['type'] == '>' and value > clause['threshold']) or
                (clause['type'] == '>=' and value >= clause['threshold']) or
                (clause['type'] == '<>' and value != clause['threshold']) or
                (clause['type'] == '=' and clause['threshold'][0] <= value <= clause['threshold'][1])
               ):
                return True

        return False


class Scenario(object):
    class NoScenarioException(Exception):
        pass

    # Map of: {'scenario ini config name': [Conditions], ...}
    condition_set = None

    @classmethod
    def _get_all_subclasses(cls):
        for subclass in cls.__subclasses__():
            yield subclass

            for subclass in subclass._get_all_subclasses():
                yield subclass

    @classmethod
    def _set_subclass_conditions(cls, config):
        for subclass in Scenario._get_all_subclasses():
            subclass.conditions = {}
            if subclass.name in config.sections():
                for item in config.items(subclass.name):
                    subclass.conditions[item[0]] = Condition(item[1])

    _subclass_conditions_set = False

    def __new__(cls, quad, scenarios_config_filename):
        if not Scenario._subclass_conditions_set:
            config = ConfigParser.ConfigParser()
            config.read(scenarios_config_filename)
            Scenario._set_subclass_conditions(config)

            Scenario._subclass_conditions_set = True

        for subclass in Scenario._get_all_subclasses():
            remaining_types = subclass.conditions.keys()
            for quad_type, quad_vaf in quad.iteritems():
                if quad_type in remaining_types and subclass.conditions[quad_type].test(quad_vaf):
                    remaining_types.remove(quad_type)
                else:
                    # Debug output:
                    #if quad_type in remaining_types:
                    #    print subclass.name, ':', quad_type, quad_vaf, 'failed', subclass.conditions[quad_type]._clauses
                    break

                if not remaining_types:
                    return object.__new__(subclass)

        raise Scenario.NoScenarioException('No matching Scenario for quad %s' % quad)
        
class GermlineMosaicAllInputs(Scenario):
    name = 'germline_mosaic_all_inputs'

class GermlineMosaicNoRnaNormal(Scenario):
	name = 'germline_mosaic_no_rna_normal'
	
class GermlineMosaicDnaOnly(Scenario):
	name = 'germline_mosaic_DNA_only'
	
class TumorInNormalAllInputs(Scenario):
	name = 'tumor_in_normal_all_inputs'
	
class TumorInNormalNoRnaNormal(Scenario):
	name = 'tumor_in_normal_no_rna_normal'

class TumorInNormalDnaOnly(Scenario):
	name = 'tumor_in_normal_DNA_only'

class RNAedAllInputs(Scenario):
    name = 'rnaed_all_inputs'

class RNAedNoRNANormal(Scenario):
    name = 'rnaed_no_rna_normal'

class TRNAedAllInputs(Scenario):
    name = 't_rnaed_all_inputs'

class VSEAllInputs(Scenario):
    name = 'vse_all_inputs'

class VSENoRNANOrmal(Scenario):
    name = 'vse_no_rna_normal'
    
class VSENormalOnly(Scenario):
    name = 'vse_normal_only'

class TVSEAllInputs(Scenario):
    name = 't_vse_all_inputs'
	
class VSETumorOnly(Scenario):
	name = 'vse_tumor_only'

class VSLAllInputs(Scenario):
    name = 'vsl_all_inputs'

class VSLNoRNANormal(Scenario):
    name = 'vsl_no_rna_normal'
    
class VSLNormalOnly(Scenario):
    name = 'vsl_normal_only'

class TVSLAllInputs(Scenario):
    name = 't_vsl_all_inputs'
	
class VSLTumorOnly(Scenario):
	name = 'vsl_tumor_only'

class LOHAltAllInputs(Scenario):
    name = 'loh_alt_all_inputs'

class LOHAltNoRNANormal(Scenario):
    name = 'loh_alt_no_rna_normal'
    
class LOHAltDNAOnly(Scenario):
    name = 'loh_alt_dna_only'
    
class LOHRefAllInputs(Scenario):
    name = 'loh_ref_all_inputs'
    
class LOHRefNoRNANormal(Scenario):
    name = 'loh_ref_no_rna_normal'
    
class LOHRefDNAOnly(Scenario):
    name = 'loh_ref_dna_only'

class SomaticAllInputs(Scenario):
    name = 'somatic_all_inputs'
    
class SomaticNoRNANormla(Scenario):
    name = 'somatic_no_rna_normal'

class SomaticDNAOnly(Scenario):
    name = 'somatic_dna_only'

class GermlineAllInputs(Scenario):
    name = 'germline_all_inputs'

class GermlineNoRnaNormal(Scenario):
    name = 'germline_no_rna_normal'

class GermlineNormalOnly(Scenario):
    name = 'germline_normal_only'

