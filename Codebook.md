# Codebook: Peace Goods for Peace

**Project Title:**  
Private Goods for Peace: Economic Provisions of Peace Agreements and the Durability of Peace  
**Code Author:** Elisa D'Amico  
**Uploaded:** 8 April 2025  

---

## Purpose  
This analysis investigates how economic provisions in peace agreements affect the duration of peace, using Cox proportional hazards models.

---

## Variables Used

### Dependent Variable
- **`duration_in_days`**: Number of days peace lasted following an agreement.

### Key Independent Variables – *Private Goods*
- **`funds_dummy_vSum22`**: Dummy for direct funds provided to ex-combatants.
- **`inclusive_dummy_vSum22`**: Dummy for inclusion/reintegration measures targeting ex-combatants.
- **`tp_all_vSum22`**: Dummy for training programs aimed at ex-combatants.
- **`PC1`**: Principal component summarizing private goods provisions.

### Key Independent Variables – *Public Goods*
- **`EpsFis`**: Dummy indicating whether fiscal federalism was included.
- **`Dev`**: Categorical variable capturing development/socio-economic reconstruction elements.
- **`NEC`**: Dummy for inclusion of a national economic plan.
- **`PC1`**: Principal component summarizing public goods provisions.

### Control Variables
- **`log_Lgt`**: Log-transformed peace agreement length.
- **`cumulative_intensity`**: Measure of total past conflict intensity.
- **`log_GDP_PPP`**: Log of GDP per capita in purchasing power parity.
- **`log_develassistance`**: Log-transformed measure of development aid received.
- **`v2x_libdem`**: V-Dem’s liberal democracy index.
- **`PrevAgmt_Bin`**: Binary variable indicating presence of previous peace agreements.
- **`post_1990`**: Dummy for post-Cold War period.
- **`post_2005`**: Dummy for post-2005 period.
- **`ddr`**: Binary indicator for disarmament, demobilization, and reintegration provisions.
- **`contype_2`** / **`contype_3`**: Conflict type indicators (categorical).
- **`year`**: Year in which agreement was signed.
- **`conflict_id`**: Unique identifier for conflict case.

---

## Notes  
- Principal components (`PC1`) are generated separately for public and private goods.
- The model is estimated using survival analysis techniques (Cox regression), accounting for clustered standard errors where applicable.
