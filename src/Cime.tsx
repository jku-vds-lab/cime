import { useState } from 'react';
import { PSEContextProvider, API, Application, PluginRegistry, createRootReducer } from 'projection-space-explorer';
import { ChemPlugin } from './Cime/ChemPlugin';
import { DatasetTabPanel } from './Overrides/DatasetTabPanel';
import { CimeAppBar } from './Overrides/CimeAppBar';
import { LineUpContext } from './LineUpContext';
import { LineUpTabPanel } from './Overrides/LineUpTabPanel';
import { AppState, cimereducers } from './State/Store';
import { Box, Slider, styled } from '@mui/material';

export const DEMO = false



PluginRegistry.getInstance().registerPlugin(new ChemPlugin())







function Cime() {
  const [context, _] = useState(new API<AppState>(null, createRootReducer(cimereducers)))

  return <PSEContextProvider context={context}>
    <Application
      config={{
        preselect: { url: 'datasets/test.sdf' }
      }}
      features={{
        disableEmbeddings: {
          tsne: true,
          forceatlas: true
        }
      }}
      overrideComponents={{
        datasetTab: DatasetTabPanel,
        appBar: CimeAppBar,
        tabs: [
          {
            name: 'lineup',
            //@ts-ignore
            tab: LineUpTabPanel,
            title: 'LineUp Integration',
            description: 'Settings for LineUp Integration',
            icon: null
          }
        ],
        detailViews: [{
          name: 'lineup',
          //@ts-ignore
          view: LineUpContext
        }]
      }} />
  </PSEContextProvider>
}

export default Cime;