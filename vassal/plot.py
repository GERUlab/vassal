""" Base class for SSA plot methods

"""

import abc
import matplotlib.pyplot as plt

class PlotSSA(object):

    __metaclass__ = abc.ABCMeta

    #---------------------------------------------------------------------------
    # Abstract methods



    # --------------------------------------------------------------------------
    # Plotting methods

    def plot(self, pltname='values', show=True, **pltkw):

        if pltname not in self._plotnames:
            names = ','.join(self._plotnames)
            raise AttributeError(
                'Unknown plot name \'{}\'. Name should be on of {}.'.format(
                    pltname, names))

        elif pltname == 'values':
            fig, ax = self._value_plot(**pltkw)

        elif pltname == 'reconstruction':
            fig, ax = self._reconstruction_plot(**pltkw)

        elif pltname == 'wcorr':
            fig, ax = self._wcorr_plot(**pltkw)

        elif pltname == 'vectors':
            fig, ax = self._vectors_plot(**pltkw)

        elif pltname == 'paired':
            fig, ax = self._paired_plot(**pltkw)

        else:
            raise NotImplementedError('{} not implemented'.format(pltname))

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])

        if show is True:
            plt.show()

        return fig, ax

    @property
    def _plotnames(self):
        names = [
            'values',
            'reconstruction',
            'wcorr',
            'vectors',
            'paired'
        ]
        return names

    def _value_plot(self, n=50, **pltkw):

        # eigenvalues
        eigenvalues = self.svd[1]  # TODO: check if needed to raise power

        #
        fig = plt.figure()
        ax = fig.gca()
        ax.semilogy(eigenvalues[:n], '-ok', markersize=4., alpha=0.5)
        ax.set_ylabel('Component Norms')
        ax.set_xlabel('Index')

        return fig, ax

    def _reconstruction_plot(self, **pltkw):

        groups = self.groups

        if len(groups) > 1:

            fig, axarr = plt.subplots(len(groups), 1, sharex=True)

            for i, g in enumerate(groups):
                ts = self._getseries(g)
                axarr[i].plot(ts, **pltkw)
                axarr[i].set_title(g)
                if i == len(groups) - 1:
                    axarr[i].set_xlabel('Time')
        else:
            pass

        return fig, axarr

    def _wcorr_plot(self, n=20, *args, **kwargs):

        wcorr = self.wcorr(components=n)

        fig = plt.figure()
        ax = fig.gca()
        im = ax.pcolor(wcorr, vmin=-1, vmax=1, cmap='PiYG')
        ax.set_aspect('equal')

        # set ticks

        ticks = np.arange(wcorr.shape[0])
        ax.set_xticks(ticks + 0.5, minor=False)
        ax.set_yticks(ticks + 0.5, minor=False)

        ax.set_xticklabels(int(x) for x in ticks)
        ax.set_yticklabels(int(x) for x in ticks)

        ax.set_title('w-correlation matrix')

        fig.colorbar(im)

        return fig, ax

    def _vectors_plot(self, n=10, **pltkw):
        """
        The rows of v are the eigenvectors of a.H a. The columns of u are the 
        eigenvectors of a a.H. For row i in v and column i in u, the 
        corresponding eigenvalue is s[i]**2.

        Parameters
        ----------
        n

        Returns
        -------

        """
        # TODO: type error
        u = self.svd[0]
        s = self.svd[1] ** 2  # TODO: check if power is needed

        fig = plt.figure()

        # grid size

        m = int(np.ceil(np.sqrt(n)))

        ax = None

        for i in range(n):
            ax = plt.subplot(m, m, i + 1, sharey=ax)
            ax.plot(u[:, i], **pltkw)
            ax.set_xticks([])
            ax.set_yticks([])

            contribution = s[i] / np.sum(s) * 100

            title = 'EV{0} ({1:.0f} %)'.format(i + 1, contribution)

            ax.set_title(title, {'fontsize': 10.})

        fig.suptitle('Eigenvectors plot')

        return fig, fig.get_axes()

    def _paired_plot(self, pairs=zip(range(0, 9), range(1, 10)), **pltkw):

        # TODO: check type pairs list of tuple of size 2

        print pairs
        u = self.svd[0]
        s = self.svd[1] ** 2  # TODO: check if power is needed

        fig = plt.figure()

        m = int(np.ceil(np.sqrt(len(pairs))))

        ax = None

        for i, j in pairs:
            ax = plt.subplot(m, m, i + 1, sharey=ax)
            ax.plot(u[:, j], u[:, i], **pltkw)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_aspect('equal')

            ssum = np.sum(s)

            contribution1 = s[i] / ssum * 100
            contribution2 = s[j] / ssum * 100

            title = 'EV{0} ({1:.0f}%) vs EV{2} ({3:.0f}%)'.format(
                i + 1,
                contribution1,
                i + 2,
                contribution2
            )

            ax.set_title(title, {'fontsize': 10.})

            fig.suptitle('Pairs of eigenvectors')

        return fig, fig.get_axes()