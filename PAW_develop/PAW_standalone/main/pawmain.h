#ifndef _PAW_MAIN_
#define _PAW_MAIN_


#ifndef CHARM_OFF
class main : public Chare {
  public:
    main(CkMigrateMessage *m) { }
    main(CkArgMsg *);
    //   ~main();
};
#else
int main(int, char *[]);
#endif

#endif
