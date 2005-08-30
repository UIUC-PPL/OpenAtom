/*-------------------------------------------------------------------------*/
/*       control_set_mol_params.c        */

 void control_set_mol_params(MDATOM_MAPS *,CPOPTS *,CPCOEFFS_INFO *,
                             CP_PARSE *,CLASS_PARSE *, MDINTRA *,
                             MDSURFACE *,FILENAME_PARSE *,FREE_PARSE *,
                             DICT_MOL *,DICT_WORD *, char *,int *,
                             int ,double ,int ,int);

/*-------------------------------------------------------------------------*/
/*       set_mol_dict.c        */

void set_molset_fun_dict(DICT_WORD *[],int *);

void set_user_base_dict(DICT_WORD *[],int *);

void set_wave_dict(DICT_WORD *[],int *, CP_PARSE *);

void set_mol_dict(DICT_WORD *[],int *,int,double,double,int);

void set_bond_free_dict(DICT_WORD *[],int *);

void set_bend_free_dict(DICT_WORD *[],int *);

void set_tors_free_dict(DICT_WORD *[],int *);

void set_rbar_free_dict(DICT_WORD *[],int *);

void set_surf_dict(DICT_WORD *[],int *);

void set_def_base_dict(DICT_WORD *[],int *);


/*-------------------------------------------------------------------------*/
/*       set_mol_params.c        */

void set_mol_params(FILENAME_PARSE*,char [], DICT_WORD [], int ,
                    CLASS_PARSE *,MDATOM_MAPS *,int [],int);

void set_wave_params(char *,char [], DICT_WORD [], int ,
                     CPOPTS *,CPCOEFFS_INFO *,CP_PARSE *);

void set_bond_free_params(char *,char [], DICT_WORD [], int ,
                          MDBOND_FREE *,FREE_PARSE *, int);

void set_bend_free_params(char *,char [], DICT_WORD [], int ,
                          MDBEND_FREE *,FREE_PARSE *, int );

void set_tors_free_params(char [],char [], DICT_WORD [], int ,
                          MDTORS_FREE *,FREE_PARSE *, int );

void set_rbar_free_params(char [],char [], DICT_WORD [], int ,
                          MDRBAR_SIG_FREE *,FREE_PARSE *, int );

void set_surf_params(char [], char [], DICT_WORD [], int ,
                          MDSURFACE *); 

void set_user_base_params(FILENAME_PARSE *,DICT_WORD [], int );

void set_def_base_params(FILENAME_PARSE *,DICT_WORD [], int );


/*-------------------------------------------------------------------------*/





