class Condition(object):
    conditions_types = ['=', '<', '<=', '>', '>=']

    def __init__(self, conditions_str):
        self.conditions = []

        for condition_str in conditions_str.split('|'):
            tokens = condition_str.strip().split(' ')
            tokens = [tok.strip() for tok in tokens]

            if len(tokens) == 2:
                threshold = float(tokens[1])
            elif len(tokens) == 3:
                threshold = (float(tokens[1]), float(tokens[2]))
            else:
                raise ValueError('Malformed conditions string - incorrect length (%s)' % conditions_str)

            if tokens[0] not in Condition.conditions_types:
                raise ValueError('Malformed conditions string - incorrect type (%s)' % conditions_str)

            self.conditions.append({'type': tokens[0], 'threshold': threshold})

    def test(self, value):
        value = float(value)
        for condition in self.conditions:
            if (
                (condition['type'] == '<' and value < condition['threshold']) or
                (condition['type'] == '<=' and value <= condition['threshold']) or
                (condition['type'] == '>' and value > condition['threshold']) or
                (condition['type'] == '>=' and value >= condition['threshold']) or
                (condition['type'] == '=' and
                    (value >= condition['threshold'][0] and value <= condition['threshold'][1]))
               ):
                return True

        return False
