Create the actual initial fields for your solid and hydraulic setup here.

At minimum, a coupled shared-mesh case usually needs user-supplied initial files for:

- displacement / solid fields from your solids4Foam setup
- `p_rgh` for the hydraulic state
- any additional fields required by your chosen solid model

This template keeps `0/` intentionally skeletal because the exact field set depends on the solid-model workflow you use.

