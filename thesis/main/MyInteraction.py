from thesis.main.BC import Integral,DirichletBC
from enum import Enum

class FieldInteractionType(Enum):
    DIRICHLET = 1
    INTEGRAL = 2
    OUTERINTERGRAL = 3

class MyInteractionTemplate:

    def __init__(self):

        self.field_quantity = None

    def get_interaction(self):
        return None

class MyFieldInteractionTemplate(MyInteractionTemplate):

    def __init__(self,field_quantity: str, field_interaction_type: FieldInteractionType):

        assert isinstance(field_quantity,str)
        assert isinstance(field_interaction_type, FieldInteractionType)

        self.field_quantity:str = field_quantity
        self.field_interaction_type: FieldInteractionType = field_interaction_type

    def get_interaction(self):

        from thesis.main.bcFunctions import cellBC as boundary_expression

        if self.field_interaction_type == FieldInteractionType.INTEGRAL:
            bc = Integral(boundary_expression, field_quantity=self.field_quantity)
            return bc


