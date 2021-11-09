import * as React from "react";
import { useState } from "react";
import {
  PSEContextProvider,
  API,
  Application,
  PluginRegistry,
  createRootReducer,
  PSEIcons,
  BaseConfig,
  FeatureConfig,
  ComponentConfig,
  DatasetType,
  setDatasetEntriesAction,
} from "projection-space-explorer";
import { ChemPlugin } from "./Cime/ChemPlugin";
import { DatasetTabPanel } from "./Overrides/DatasetTabPanel";
import { CimeAppBar } from "./Overrides/CimeAppBar";
import { LineUpContext } from "./LineUpContext";
import { LineUpTabPanel } from "./Overrides/LineUpTabPanel";
import { AppState, CIMEReducers } from "./State/Store";
import merge from 'lodash/merge'
import clone from "lodash/clone";

export const DEMO = false;

PluginRegistry.getInstance().registerPlugin(new ChemPlugin());


export type CIMEAppProps = {
  config?: BaseConfig
  features?: FeatureConfig
  overrideComponents?: ComponentConfig
}

export const DEFAULT_CIME_APP_CONFIG = {
  config: {
    preselect: {
      initOnMount: false
    }
  },
  features: {
    disableEmbeddings: {
      tsne: true,
      forceatlas: true,
    },
  },
  overrideComponents: {
    datasetTab: DatasetTabPanel,
    appBar: CimeAppBar,
    tabs: [
      {
        name: "lineup",
        tab: LineUpTabPanel,
        title: "LineUp Integration",
        description: "Settings for LineUp Integration",
        icon: PSEIcons.PseLineup,
      },
    ],
    detailViews: [
      {
        name: "lineup",
        view: LineUpContext,
      },
    ],
  }
}


fetch('/jku-vds-lab/cime/bundle.js', {
  mode: 'no-cors'
}).then((response) => response.text())
.then(json => {
  console.log(json)
})
.catch(() => {
  console.log("wrong")
})


export function CIMEApp(props: CIMEAppProps) {
  const api = new API<AppState>(null, createRootReducer(CIMEReducers))

  api.store.dispatch(setDatasetEntriesAction([
    {
      display: "Chess: 190 Games",
      path: "domain_5000_experiments_0.csv",
      type: DatasetType.None
    },
    {
      display: "Chess: 450 Games",
      path: "cube1x2.csv",
      type: DatasetType.Rubik
    }
  ]))
  
  const [context] = useState(
    api
  );

  const [merged, setMerged] = React.useState(merge(clone(DEFAULT_CIME_APP_CONFIG), props))

  React.useEffect(() => {
    const cimeAppConfig = clone(DEFAULT_CIME_APP_CONFIG)

    setMerged(merge(cimeAppConfig, props))
  }, [props.config, props.features, props.overrideComponents])

  return (
    <PSEContextProvider context={context}>
      <Application
      config={merged.config}
      features={merged.features}
      //@ts-ignore
      overrideComponents={merged.overrideComponents}
      />
    </PSEContextProvider>
  );
}