import ConfigParser
from collections import defaultdict
from maf_types import MafTypes

SCENARIOS_CONFIG = 'scenarios.ini'


class Event(object):
    germline_mosaic = 'germline_mosaic'
    tumor_in_normal = 'tumor_in_normal'
    vse_all_inputs = 'vse_all_inputs'
    t_vse = 't_vse'
    vse = 'vse'
    vsl = 'vsl'
    t_vsl = 't_vsl'
    loh_alt = 'loh_alt'
    loh_ref = 'loh_ref'
    germline = 'germline'
    rnaed = 'rnaed'
    t_rnaed = 't_rnaed'
    somatic = 'somatic'
    unclassified = 'unclassified'

    @classmethod
    def all(cls):
        return [Event.germline_mosaic,
                Event.tumor_in_normal,
                Event.vse_all_inputs,
                Event.vse,
                Event.t_vse,
                Event.vsl,
                Event.t_vsl,
                Event.loh_alt,
                Event.loh_ref,
                Event.germline,
                Event.rnaed,
                Event.t_rnaed,
                Event.somatic,
                Event.unclassified]


class Condition(object):
    comparators = ['=', '<', '<=', '>', '>=', '<>']

    def __init__(self, condition_str):
        self._clauses = []

        for clause in condition_str.split('|'):
            tokens = clause.strip().split(' ')
            tokens = [tok.strip() for tok in tokens]
            comparator = tokens[0]

            if len(tokens) == 2:
                thresholds = [float(tokens[1])]
                if comparator == '=':
                    raise ValueError(
                        'Equality must be specified using two thresholds (= 0.0 0.0, for example, instead of = 0.0): {}'
                            .format(condition_str))
            elif len(tokens) == 3:
                thresholds = [float(tokens[1]), float(tokens[2])]
            else:
                raise ValueError('Malformed conditions string - incorrect length ({})'.format(condition_str))

            if comparator not in Condition.comparators:
                raise ValueError('Malformed conditions string - incorrect comparator ({})'.format(condition_str))

            for t in thresholds:
                if t > 1:
                    raise ValueError('Threshold {} is larger than 1'.format(t))
                if t < 0:
                    raise ValueError('Threshold {} is less than 0'.format(t))
            self._clauses.append({'type': tokens[0],
                                  'threshold': thresholds})

    def test(self, value):
        value = float(value)
        for clause in self._clauses:
            if (
                (clause['type'] == '<' and value < clause['threshold'][0]) or
                (clause['type'] == '<=' and value <= clause['threshold'][0]) or
                (clause['type'] == '>' and value > clause['threshold'][0]) or
                (clause['type'] == '>=' and value >= clause['threshold'][0]) or
                (clause['type'] == '<>' and value != clause['threshold'][0]) or
                (clause['type'] == '=' and clause['threshold'][0] <= value <= clause['threshold'][1])
               ):
                return True

        return False


class ScenarioCalculator(object):
    class NoScenarioException(Exception):
        pass

    class UnrecognizedEventException(Exception):
        pass

    analysis_decision_tree = defaultdict(lambda: defaultdict(list))

    def __build_decision_tree(self, scenarios_config_filename):
        """Config sections should be labeled <analysis_type>.<gl_vaf_status>.<priority>.<event>. Figure out what the
        set of possible analysis types is based on the config sections present in the .ini file, and build the decision
        tree for triaging requests evaluate VAF groupings."""
        config = ConfigParser.ConfigParser()
        config.read(scenarios_config_filename)

        for c in config.sections():
            analysis_type, germline_vaf, priority, event = c.split('.')
            if event not in Event.all():
                raise self.UnrecognizedEventException("Event {} not expected".format(event))
            priority = int(priority)

            conditions = {}

            for maf_type in MafTypes.all():
                try:
                    condition_for_type = config.get(c, maf_type)
                    conditions[maf_type] = Condition(condition_for_type)
                except ConfigParser.NoOptionError:
                    continue

            self.analysis_decision_tree[analysis_type][germline_vaf].append({
                'priority': priority,
                'event': event,
                'conditions': conditions
            })

    def __init__(self, scenarios_config_filename):
        self.__build_decision_tree(scenarios_config_filename)

    def __get_possible_categories(self, analysis_type, vaf_values):
        # Navigate decision tree first by analysis type
        categories_for_analysis_type = self.analysis_decision_tree.get(analysis_type, None)
        if not categories_for_analysis_type:
            raise self.NoScenarioException("No scenario found for analysis type {}".format(analysis_type))

        # Further navigate tree by germline VAF status
        if 'dna_normal' in vaf_values.keys():
            germline_vaf = vaf_values['dna_normal']
            if germline_vaf > 0:
                return categories_for_analysis_type['gl>0']
            elif germline_vaf == 0:
                return categories_for_analysis_type['gl=0']
        else:
            return categories_for_analysis_type['no_gl']

    def categorize(self, analysis_type, vaf_values):
        possible_categories = self.__get_possible_categories(analysis_type, vaf_values)
        # For each category, by priority order, check whether the VAFs match for the event. Return the first event
        # for which the VAFs fit the criteria.
        for category in sorted(possible_categories, key=lambda x: x['priority']):
            category_fits = True
            for maf_type, condition in category.get('conditions').items():
                if maf_type in vaf_values:
                    # If any condition fails for this category, we will not select this category
                    if not condition.test(vaf_values.get(maf_type)):
                        category_fits = False
            # If all conditions have passed for this category, we will select this category
            if category_fits:
                return category.get('event')

        raise self.NoScenarioException("No scenario found for analysis type {} and vaf values {}".format(analysis_type,
                                                                                                         vaf_values))
