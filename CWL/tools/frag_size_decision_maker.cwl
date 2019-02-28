cwlVersion: v1.0
class: ExpressionTool
requirements:
  InlineJavascriptRequirement: {}
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 1000

inputs:
  user_def_fragment_size:
    type: int?
  cc_fragment_size:
    type: int?

expression: |
  ${
      var fragment_size = inputs.user_def_fragment_size;
      if( fragment_size == null ){
          fragment_size = cc_fragment_size;
      }
      return { "fragment_size": fragment_size }
  }

outputs:
  fragment_size:
    type: int?