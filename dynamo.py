
# Imports:
from traits.api import HasPrivateTraits, Float, Range, Int

from traitsui.api import View, Group, Item, Label, RangeEditor


class DynamicRangeEditor(HasPrivateTraits):
    """ Defines an editor for dynamic ranges (i.e. ranges whose bounds can be
        changed at run time).
    """

    # The value with the dynamic range:
    value = Float

    # This determines the low end of the range:
    low = Range(0.0, 10.0, 0.0)

    # This determines the high end of the range:
    high = Range(20.0, 100.0, 20.0)

    # An integer value:
    int_value = Int

    # This determines the low end of the integer range:
    int_low = Range(0, 10, 0)

    # This determines the high end of the range:
    int_high = Range(20, 100, 20)

    # Traits view definitions:
    traits_view = View(

        # Dynamic simple slider demo:
        Group(
            Item('value',
                 editor=RangeEditor(low_name='low',
                                    high_name='high',
                                    format='%.1f',
                                    label_width=28,
                                    mode='auto')
                 ),
            '_',
            Item('low'),
            Item('high'),
            '_',
            Label('Move the Low and High sliders to change the range of '
                  'Value.'),
            label='Simple Slider'
        ),

        # Dynamic large range slider demo:
        Group(
            Item('value',
                 editor=RangeEditor(low_name='low',
                                    high_name='high',
                                    format='%.1f',
                                    label_width=28,
                                    mode='xslider')
                 ),
            '_',
            Item('low'),
            Item('high'),
            '_',
            Label('Move the Low and High sliders to change the range of '
                  'Value.'),
            label='Large Range Slider'
        ),

        # Dynamic spinner demo:
        Group(
            Item('int_value',
                 editor=RangeEditor(low=0,
                                    high=20,
                                    low_name='int_low',
                                    high_name='int_high',
                                    format='%d',
                                    is_float=False,
                                    label_width=28,
                                    mode='spinner')
                 ),
            '_',
            Item('int_low'),
            Item('int_high'),
            '_',
            Label('Move the Low and High sliders to change the range of '
                  'Value.'),
            label='Spinner'
        ),
        title='Dynamic Range Editor Demonstration',
        buttons=['OK'],
        resizable=True
    )


# Create the demo:
demo = DynamicRangeEditor()

# Run the demo (if invoked from the command line):
if __name__ == '__main__':
    demo.configure_traits()