import numbers, os

from pyaeroopt.interface import CodeInterface
from pyaeroopt.io        import InputBlock, InputFile
from pyaeroopt.util.hpc_util import execute_str, execute_code
from pyaeroopt.util.sdesign_util import tmp_fname, tmp_femesh
from pyaeroopt.util.sdesign_util import prepare_directory, clean_directory

class Sdesign(CodeInterface):
    """
    An object to facilitate interfacing to SDESIGN.
    """
    def __init__(self, **kwargs):

        # Constructor of base class, check fields, and extract anticipated input
        super(Sdesign, self).__init__(**kwargs)

        if self.bin is None: self.bin = os.path.expandvars('$SDESIGN')


class SdesignInputFile(InputFile):
    def __init__(self, fname, blocks, log='sdesign.tmp.log', **kwargs):
        super(SdesignInputFile, self).__init__(fname, blocks)
        self.sep = SdesignInputBlock.line_break
        self.log = log
        for kwarg in kwargs:
            setattr(self, kwarg, kwargs[kwarg])

    def write(self):
        InputFile.write(self)
        f=open(self.fname, 'a+')
        f.write('\nEND')
        f.close()

class SdesignInputBlock(InputBlock):
    def __init__(self, name, *args):
        super(SdesignInputBlock, self).__init__(name, *args)

    def write(self, fname, indent_level):

        blockvars = vars(self)
        indent = ''
        f=open(fname,'a+')
        f.write(indent*indent_level+self.name)

        # Formatting for FEMESH vs. rest of blocks
        if self.name == 'FEMESH': f.write('  ')
        else:                     f.write('\n')

        # Print all of the properties of this block
        for prop in self.props:
            if prop == 'name': continue

            # Allow comments as property blocks of the form ['# ....'] in any property block
            if prop[0] == '#':
                f.write(prop + '\n')
                continue

            # DEFINE block properties should have the form ['x1', 'expression for x1']
            if self.name == 'DEFINE':
                f.write(indent*(indent_level+1)+prop+' = '
                                                     +str(blockvars[prop])+'\n')
                continue
            if self.name == 'FEMESH':
                f.write(indent*(indent_level+1)+'"'+str(prop)+'"\n')
                continue   
            if isinstance(blockvars[prop],(str,numbers.Number)):
                f.write(indent*(indent_level+1)+str(prop + '    '+blockvars[prop] + '\n'))
            if isinstance(blockvars[prop],(list,tuple)):
                f.write(indent*(indent_level+1)+str(prop) + '    ')
                for x in blockvars[prop]:
                    if isinstance(x,(str,numbers.Number)):
                        f.write(indent*(indent_level+1)+str(x)+'    ')
                    if isinstance(x,(list,tuple)):
                        for y in x:
                            if isinstance(y,(str,numbers.Number)):
                                f.write(indent*(indent_level+1)+str(y)+'    ')
                        f.write('\n')
                f.write('\n')

        f.close()

    @staticmethod
    def line_break(fname):
        f=open(fname,'a'); f.write('\n'); f.close();
