Create the actual initial fields for your solid and hydraulic setup here.

For a variably saturated shared-mesh case, you will usually need user-supplied initial files for:

- displacement / solid fields from your solids4Foam setup
- either `p_rgh` or `pHead`, depending on the hydraulic formulation
- any additional solid fields required by your chosen constitutive law

The fields `S` and `n` are normally handled by the coupled poromechanics stack once the case dictionaries are consistent.

