
/* $Id$ */

#ifndef READTNC2_H
#define READTNC2_H

#ifdef __cplusplus
extern "C" {
#endif
  
  void create_read_tnc2(const char *file_name);
  void delete_read_tnc2();
  
  int snapshot_size();
  int n_atoms();
  
  int n_active_atoms();
  int n_total_qm_atoms();
  
  const int *atom_types();
  const char *atom_names();
  const double *snapshot(const int i);
  
  const int *active_atom_list();
  const double *qm_esp_charges();
  const double *gradients();
  const double *velocities();
  
  int qm_esp_charges_saved();
  int gradients_saved();
  int velocities_saved();
  
#ifdef __cplusplus
}
#endif

#endif /* READTNC2_H */
