'''
@author: Yuxiang Chen
'''

import signal
import sys

class MPI:
    """docstring for MPI"""
    def __init__(self):
        try: # check if the library is installed or not
            from mpi4py import MPI
        except:
            raise Exception("mpi4py library is not installed!")

        self.comm = MPI.COMM_WORLD
        self.size = self.comm.Get_size()
        self.rank = self.comm.Get_rank()

        self._begun = False

        #User termination on Ctrl-C will be caught and send to workers.
        signal.signal(signal.SIGINT, self.signal_handler)
        signal.signal(signal.SIGTERM, self.signal_handler)
        #signal.signal(signal.SIGKILL, self.signal_handler)

    def signal_handler(self, sig, frame):
        self.comm.Abort()
        print('Sent abort signal to comm_world')
        sys.exit(0)

    def is_master(self):
        return self.rank == 0
    

    def _split_seq(self, seq, size):
        if len(seq) == size:
            return seq
        elif len(seq) > size:
            new_seq = [None] * size
            n, N = 0, len(seq)
            for i in range(size):
                l = N // size + (N % size > i)
                new_seq[i] = seq[n:n+l]
                n += l
            return new_seq
        else:
            new_seq = seq + [None] * (size-len(seq))
            return new_seq


    def _merge_seq(self, seq, size):
        if len(seq) == size:
            return seq
        elif len(seq) > size:
            return seq[:size]
        else:
            new_seq = []
            for s in seq:
                new_seq += s
            assert len(new_seq) == size
            return new_seq


    def begin(self):
        """For master and worker.
        """
        if not self._begun:
            self._begun = True
            if self.is_master():
                return
            else:
                self.standby()


    def end(self):
        """For master only.
        """
        if self._begun:
            if self.is_master():
                msg = [None] * self.size
                self.comm.scatter(msg, root=0)
        self._begun = False


    def standby(self):
        """For worker only.
        """
        if not self._begun:
            raise Exception("MPI is not begun!")

        if self.is_master():
            return

        while True:
            try:
                # get the msg from master
                msg = self.comm.scatter(None, root=0)
                if msg is None: # master send end msg
                    break
                # the first part is func and second part is data
                assert len(msg) == 2
                func = msg[0]
                data = msg[1]
                if data is None: # no job to do
                    res = None
                elif data.__class__ == list:
                    res = []
                    for d in data:
                        if d.__class__ == tuple:
                            res.append(func(*d))
                        else:
                            res.append(func(d))
                else:
                    if data.__class__ == tuple:
                        res = func(*data)
                    else:
                        res = func(data)

                # send back the result
                self.comm.gather(res, root=0)

            except Exception as e:
                print(e)
                self.comm.Abort()

        # get end msg, terminate
        import sys
        sys.exit()


    def parfor(self, func, data, verbose=False):
        """For master only.
        """
        if not self._begun:
            raise Exception("MPI has not been initialized!")

        if not self.is_master():
            return

        assert func
        assert data.__class__ == list
        assert len(data)

        ddata = self._split_seq(data, self.size)
        msg = [(func, d) for d in ddata]

        try:
            sdata = self.comm.scatter(msg, root=0) # for generic Python objects

            # the master also has to do the job :(
            dd = sdata[1]
            if dd is None:
                res = None
            elif dd.__class__ == list:
                res = []
                for d in dd:
                    if d.__class__ == tuple: # multi args
                        res.append(func(*d))
                    else:
                        res.append(func(d))
            else:
                if dd.__class__ == tuple:
                    res = func(*dd)
                else:
                    res = func(dd)

            # gather results
            all_res = self.comm.gather(res, root=0)
            all_res = self._merge_seq(all_res, len(data))

        except Exception as e:
            print(e)
            self.comm.Abort()

        return all_res


