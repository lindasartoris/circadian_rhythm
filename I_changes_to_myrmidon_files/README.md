### 1. Average ant length for each set-up 
Manual orientation and body length for one replicate per tracking system in fort studio

These are the target files for step 2

### 2. Capsule assignment for manually oriented files
One capsule for all interactions (interaction network) and one for grooming

Based on Adriano’s capsules:

_Interactions_: “CapsuleDef2018”
 
_Grooming_: “CapsuleDef3”
 
Run **“1_CAPSULE_ASSIGNEMENT_Linda.R”** once for interaction, once for grooming

Extract the average ant sizes 

“Mean_ant_length_colonies.txt”

This gives us two myrmidon files per manually oriented target file 

“[…]_manual_orientation_CapsuleDef2018.myrmidon”

“[…]_manual_orientation_CapsuleDef3.myrmidon”

### 3. Remove corrupted files
Because of a power cut & a problem with starting the tracking (which caused issues when producing the auto_orient files) some folders or files were corrupted and needed to be removed.

### 4. Create the automated ants with interaction capsule
First, only do this for the interaction capsule (Def2018) because we still need to create metadata (assign dead ants, exposed ants etc.). Metadata will be created on the Def2018.myrmidon files. Next, we remove the capsules and create new files for Def3 that will already have all the metadata.

Based on the “[…]_manual_orientation_CapsuleDef2018.myrmidon” files that exist for one replicate per tracking system

Run **“2_AUTO_ORIENT_Linda.R”**

 “[…]_automatically_oriented_CapsuleDef2018.myrmidon”

### 5. Add metadata
#### 5.1 create myrmidon files for forager_ID tracking
        
Move all the extra forager tracking files (those for lost tags, re-tags etc) to a folder “forager_IDs_extra”

-> left with one tracking file per colony

Adapt the script from “2_AUTO_ORIENT_Linda.R” to create base myrmidon files in “vital_Linda.R”

In **“vital_Linda.R”** run the section 5.1

“foragerID_c[xx].myrmidon”

#### 5.2.1 create meta data keys
Adapt the script from “vital_Daniel.R” to create the meta data keys (is.queen, is.alive, exposed, …) for the “automatically_oriented_Def2018.myrmidon” files

In **“vital_Linda.R”** run the section 5.2 (5.2.1 meta data keys are created)

#### 5.2.2 add treated forager IDs
This is included in “vital_Linda.R” section 5.2 (5.2.2 exposed forager information is added)

In that same loop the meta data key is changed to “exposed” == TRUE for the exposed foragers

Each colony has a separate tracking folder of the treated foragers with the myrmidon file (“foragerID_c[xx].myrmidon”) that can be matched by their colony ID (cXX)

“[…]_automatically_oriented_CapsuleDef2018_metadata.myrmidon”

#### 5.3 manually in myrmidon files: dead ants & queen

Ants that either died during the experiment, disappeared or lost their tag will be excluded from the analysis because the same ants have to be present in the pre- & post networks

Open the last video in the visualisation workspace and check if ant is still moving, else track her last appearance (statistics workspace for time last seen) and change meta data “still.alive” to false

Save as:

“[…]_automatically_oriented_CapsuleDef2018_metadata_deaths.myrmidon”

After doing this, include step 6. for the queen (so myrmidon files only have to be opened once)

### 6. Queen size and capsule
Manually orient the queen by setting her body length and draw the capsules (body & head)

“[…]_automatically_oriented_CapsuleDef2018_metadata_deaths_q.myrmidon”

### 7. Record re-tagging/re-orientation of ants
In myrmidon files (from notes) – only ever re-tagged with the same tag so only re-orientation was needed

E.g. use the time point when the forager was removed from the arena as end point (recorded in tag_loss_time.txt)

Also check auto orientation of other ants bc some replicates seem to have had problems (prideaux_c03)

“[…]_automatically_oriented_CapsuleDef2018_metadata_deaths_q.myrmidon”

### 8. Zoning
Create the spaces for nest & foraging

Use script **“3_CLONE_ZONES_Linda.R”** (adapted from the one for Beki by Nathalie)

Create the zones once manually for replicate alleline_c05 in workspace “zoning” and use this as a template “[…]_zones_Linda_2023.myrmidon”

Then clone (with script) and adjust manually in fort-studio where needed

“[…]_automatically_oriented_CapsuleDef2018_metadata_deaths_q_zones.myrmidon”

### 9. Create the grooming capsules
Do this once all the relevant meta data and other information is added to the myrmidon files so nothing has to be repeated for the second capsule files:

Caution: we don’t want to re-create the ants (because of the metadata) nor to auto_orient them again (because of some re-tagging)! Only remove the Def2018 capsules and create Def3 grooming capsules

Use script **“4_CLONING_CAPSULES_Linda.R”** (adapted from Luke’s script “CAPSULE_REASSIGNMENT_2811.R”)

“[…]_automatically_oriented_CapsuleDef3_metadata_deaths_q_zones.myrmidon”
