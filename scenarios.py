import ConfigParser

SCENARIOS_CONFIG = 'scenarios.ini'

class Condition(object):
    conditions_types = ['=', '<', '<=', '>', '>=']

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


class RNAedAllInputs(Scenario):
    name = 'rnaed_all_inputs'


class RNAedNormalOnly(Scenario):
    name = 'rnaed_normal_only'


class TRNAedAllInputs(Scenario):
    name = 't_rnaed_all_inputs'


class TRNAedNoDNANormal(Scenario):
    name = 't_rnaed_no_dna_normal'


class VSEAllInputs(Scenario):
    name = 'vse_all_inputs'


class VSENormalOnly(Scenario):
    name = 'vse_normal_only'


class TVSEAllInputs(Scenario):
    name = 't_vse_all_inputs'


class TVSENoDNANormal(Scenario):
    name = 't_vse_no_dna_normal'


class VSLAllInputs(Scenario):
    name = 'vsl_all_inputs'


class VSLNormalOnly(Scenario):
    name = 'vsl_normal_only'


class TVSLAllInputs(Scenario):
    name = 't_vsl_all_inputs'


class TVSLNoDNANormal(Scenario):
    name = 't_vsl_no_dna_normal'


class LOHAltAllInputs(Scenario):
    name = 'loh_alt_all_inputs'


class LOHAltDNAOnly(Scenario):
    name = 'loh_alt_dna_only'


class LOHCNAllInputs(Scenario):
    name = 'loh_cn_all_inputs'


class LOHCNDNAOnly(Scenario):
    name = 'loh_cn_dna_only'


class SomaticAllInputs(Scenario):
    name = 'somatic_all_inputs'


class SomaticDNAOnly(Scenario):
    name = 'somatic_dna_only'
