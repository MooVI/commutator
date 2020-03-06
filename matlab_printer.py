from sympy.printing.str import StrPrinter
from sympy.printing.precedence import precedence, PRECEDENCE
from sympy.core import S
from sympy.printing.precedence import precedence, PRECEDENCE
import re

def multiple_replace(string, rep_dict):
    pattern = re.compile("|".join([re.escape(k) for k in rep_dict.keys()]), re.M)
    return pattern.sub(lambda x: rep_dict[x.group(0)], string)

class MatlabStrPrinter(StrPrinter):
    """
    Examples of how to customize the StrPrinter for both a SymPy class and a
    user defined class subclassed from the SymPy Basic class.
    """

    def _print_PolyElement(self, poly):
        return poly.str(self, PRECEDENCE, "%s.^%s", "*")

    def _print_Poly(self, expr):
        terms, gens = [], [ self._print(s) for s in expr.gens ]

        for monom, coeff in expr.terms():
            s_monom = []

            for i, exp in enumerate(monom):
                if exp > 0:
                    if exp == 1:
                        s_monom.append(gens[i])
                    else:
                        s_monom.append(gens[i] + ".^%d" % exp)

            s_monom = "*".join(s_monom)

            if coeff.is_Add:
                if s_monom:
                    s_coeff = "(" + self._print(coeff) + ")"
                else:
                    s_coeff = self._print(coeff)
            else:
                if s_monom:
                    if coeff is S.One:
                        terms.extend(['+', s_monom])
                        continue

                    if coeff is S.NegativeOne:
                        terms.extend(['-', s_monom])
                        continue

                s_coeff = self._print(coeff)

            if not s_monom:
                s_term = s_coeff
            else:
                s_term = s_coeff + "*" + s_monom

            if s_term.startswith('-'):
                terms.extend(['-', s_term[1:]])
            else:
                terms.extend(['+', s_term])

        if terms[0] in ['-', '+']:
            modifier = terms.pop(0)

            if modifier == '-':
                terms[0] = '-' + terms[0]

        format = expr.__class__.__name__ + "(%s, %s"

        from sympy.polys.polyerrors import PolynomialError

        try:
            format += ", modulus=%s" % expr.get_modulus()
        except PolynomialError:
            format += ", domain='%s'" % expr.get_domain()

        format += ")"

        return format % (' '.join(terms), ', '.join(gens))

    def _print_Pow(self, expr, rational=False):
        PREC = precedence(expr)

        #if expr.base is RootOfUnity:
        #    return 'Exp['+str(2*expr.base.n)+ 'I Pi /'+ str(expr.base.n) + ']'

        if expr.exp is S.Half and not rational:
            return "sqrt(%s)" % self._print(expr.base)

        if expr.is_commutative:
            if -expr.exp is S.Half and not rational:
                # Note: Don't test "expr.exp == -S.Half" here, because that will
                # match -0.5, which we don't want.
                return "1/sqrt(%s)" % self._print(expr.base)
            if expr.exp is -S.One:
                # Similarly to the S.Half case, don't test with "==" here.
                return '1/%s' % self.parenthesize(expr.base, PREC)

        e = self.parenthesize(expr.exp, PREC)
        if self.printmethod == '_sympyrepr' and expr.exp.is_Rational and expr.exp.q != 1:
            # the parenthesized exp should be '(Rational(a, b))' so strip parens,
            # but just check to be sure.
            if e.startswith('(Rational'):
                return '%s.^%s' % (self.parenthesize(expr.base, PREC), e[1:-1])
        return '%s.^%s' % (self.parenthesize(expr.base, PREC), e)

    def _print_MatPow(self, expr):
        PREC = precedence(expr)
        return '%s^%s' % (self.parenthesize(expr.base, PREC),
                         self.parenthesize(expr.exp, PREC))

    def _print_RootOfUnity(self, expr):
        #return 'Exp[2 I Pi /'+ str(expr.n) + ']'
        return "sym('w')"

    def _print_list(self, expr):
        return "[%s]" % self.stringify(expr, " ")

    def _print_MatrixBase(self, expr):
        if expr.rows == 0 or expr.cols == 0:
            return '[]'
        l = expr.tolist()
        return "[" +";".join(self._print(sl) for sl in l)+"]"

    _print_SparseMatrix = \
        _print_MutableSparseMatrix = \
        _print_ImmutableSparseMatrix = \
        _print_Matrix = \
        _print_DenseMatrix = \
        _print_MutableDenseMatrix = \
        _print_ImmutableMatrix = \
        _print_ImmutableDenseMatrix = \
        _print_MatrixBase

    def _print_MatrixElement(self, expr):
        return self._print(expr.parent) + '[%s %s]'%(expr.i, expr.j)

    def _print_MatrixSlice(self, expr):
        def strslice(x):
            x = list(x)
            if x[2] == 1:
                del x[2]
            if x[1] == x[0] + 1:
                del x[1]
            if x[0] == 0:
                x[0] = ''
            return ':'.join(map(self._print, x))
        return (self._print(expr.parent) + '[' +
                strslice(expr.rowslice) + ' ' +
                strslice(expr.colslice) + ']')

    def _print_Pi(self, expr):
        return 'pi'

    def _print_ImaginaryUnit(self, expr):
        return 'i'

    def _print_Integer(self, expr):
        return "sym('"+str(expr.p)+"')"

    def _print_int(self, expr):
        return "sym('"+str(expr)+"')"

    def _print_mpz(self, expr):
        return "sym('"+str(expr)+"')"

    def _print_Rational(self, expr):
        if expr.q == 1:
            return "sym('"+str(expr.p)+"')"
        else:
            return "sym('%s')/sym('%s')" % (expr.p, expr.q)

    def _print_PythonRational(self, expr):
        if expr.q == 1:
            return "sym('"+str(expr.p)+"')"
        else:
            return "%d/%d" % (expr.p, expr.q)

    def _print_Fraction(self, expr):
        if expr.denominator == 1:
            return "sym('"+str(expr.numerator)+"')"
        else:
            return "sym('%s')/sym('%s')" % (expr.numerator, expr.denominator)


matlabprinter = MatlabStrPrinter()

def matlabstr(expr):
    """Returns the expression as a string.

    For large expressions where speed is a concern, use the setting
    order='none'.

    Examples
    ========

    >>> from sympy import symbols, Eq, sstr
    >>> a, b = symbols('a b')
    >>> sstr(Eq(a + b, 0))
    'Eq(a + b, 0)'
    """
    return matlabprinter.doprint(expr)
