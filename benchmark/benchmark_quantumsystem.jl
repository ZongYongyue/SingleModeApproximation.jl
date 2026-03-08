using Pkg
Pkg.activate(".")
using QuantumLattices: Lattice as QLattice, bonds as qbonds, Hilbert, expand, Fock, OperatorGenerator, Hopping, Hubbard,InterOrbitalInterSpin,InterOrbitalIntraSpin,SpinFlip,PairHopping
using MeanFieldTheories: Lattice as SLattice, Dof, SystemDofs, QN, bonds as sbonds, c, cdag, generate_onebody, generate_twobody

#QuantumLattices#############################################
lattice = QLattice([0.0], [1.0]);

# define the internal degrees of freedom, i.e., the single-orbital spin-1/2 fermionic algebra
hilbert = Hilbert(site=>Fock{:f}(2, 2) for site=1:length(lattice))

# define the terms
t = Hopping(:t, -1.0, 1)
U = Hubbard(:U, 6.0)
Up = InterOrbitalInterSpin(:Up, 3.6)
J = 1.2
UpJ = InterOrbitalIntraSpin(:UpJ, 3.6-J)
J1 = SpinFlip(:J1, J)
J2 = PairHopping(:J2, J)

us = expand(OperatorGenerator(qbonds(lattice, 1), hilbert, (U,)))
ups = expand(OperatorGenerator(qbonds(lattice, 1), hilbert, (Up,)))
ujs = expand(OperatorGenerator(qbonds(lattice, 1), hilbert, (UpJ,)))
js1 = expand(OperatorGenerator(qbonds(lattice, 1), hilbert, (J1,)))
js2 = expand(OperatorGenerator(qbonds(lattice, 1), hilbert, (J2,)))



#MeanFieldTheories#############################################
dofs = SystemDofs([
    Dof(:site, 2),
    Dof(:orbital, 2),
    Dof(:spin, 2, [:up, :down ])])
unit = SLattice(
    [Dof(:site, 1), ],
    [QN(site=1), ],
    [[0.0, 0.0], ];
)
latt = SLattice(unit, [[1.0, 0.0], [0.0, 1.0]], (2, 1))

sus = generate_twobody(dofs, sbonds(latt, (:o,:o), 0), 
    (delta, qn1, qn2, qn3, qn4) ->
        (qn1.orbital == qn2.orbital == qn3.orbital == qn4.orbital) &&
        (qn1.spin, qn2.spin, qn3.spin, qn4.spin) == (1, 1, 2, 2) ? 6.0 : 0.0
    ,
    order = (cdag, :i, c, :i, cdag, :i, c, :i)
    )
sups = generate_twobody(dofs, sbonds(latt, (:o,:o), 0), 
    (delta, qn1, qn2, qn3, qn4) ->
        (qn1.orbital, qn3.orbital) == (qn2.orbital, qn4.orbital) &&
        (qn1.orbital < qn3.orbital) &&
        (qn1.spin, qn3.spin) == (qn2.spin, qn4.spin) &&
        (qn1.spin !== qn3.spin) ? 3.6 : 0.0
    ,
    order = (cdag, :i, c, :i, cdag, :i, c, :i)
    )
sujs = generate_twobody(dofs, sbonds(latt, (:o,:o), 0), 
    (delta, qn1, qn2, qn3, qn4) ->
        (qn1.orbital, qn3.orbital) == (qn2.orbital, qn4.orbital) &&
        (qn1.orbital < qn3.orbital) &&
        (qn1.spin, qn3.spin) == (qn2.spin, qn4.spin) &&
        (qn1.spin == qn3.spin) ? 2.4 : 0.0
    ,
    order = (cdag, :i, c, :i, cdag, :i, c, :i)
    )
sjs1 = generate_twobody(dofs, sbonds(latt, (:o,:o), 0), 
    (delta, qn1, qn2, qn3, qn4) ->
        (qn1.orbital, qn2.orbital) == (qn3.orbital, qn4.orbital) &&
        (qn1.orbital !== qn2.orbital) &&
        (qn1.spin, qn2.spin, qn3.spin, qn4.spin) == (1,2,2,1) ? 1.2 : 0.0
    ,
    order = (cdag, :i, cdag, :i, c, :i, c, :i)
    )
sjs2 = generate_twobody(dofs, sbonds(latt, (:o,:o), 0), 
    (delta, qn1, qn2, qn3, qn4) ->
        (qn1.orbital, qn3.orbital) == (qn2.orbital, qn4.orbital) &&
        (qn1.orbital !== qn3.orbital) &&
        (qn1.spin, qn2.spin, qn3.spin, qn4.spin) == (1,2,2,1) ? 1.2 : 0.0
    ,
    order = (cdag, :i, cdag, :i, c, :i, c, :i)
    )