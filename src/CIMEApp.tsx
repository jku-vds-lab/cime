import * as React from "react";
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
} from "projection-space-explorer";
import { ChemPlugin } from "./Cime/ChemPlugin";
import { DatasetTabPanel } from "./Overrides/DatasetTabPanel";
import { CimeAppBar } from "./Overrides/CimeAppBar";
import { LineUpContext } from "./LineUpContext";
import { LineUpTabPanel } from "./Overrides/LineUpTabPanel";
import { AppState, CIMEReducers } from "./State/Store";
import merge from 'lodash/merge'
import clone from "lodash/clone";
import cloneDeep from 'lodash/cloneDeep'
import { ConnectedComponent } from "react-redux";

export const DEMO = false;

PluginRegistry.getInstance().registerPlugin(new ChemPlugin());


export type CIMEAppProps = {
  config?: BaseConfig;
  features?: FeatureConfig;
  overrideComponents?: ComponentConfig;
  pseRef?: any;
  providePSEContext?: boolean;
}

export const DEFAULT_CIME_APP_CONFIG: CIMEAppProps = {
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
        // @ts-ignore TODO: @moritz error after I added the correct typing for DEFAULT_CIME_APP_CONFIG
        tab: LineUpTabPanel,
        title: "LineUp Integration",
        description: "Settings for LineUp Integration",
        icon: PSEIcons.PseLineup,
      },
    ],
    detailViews: [
      {
        name: "lineup",
        // @ts-ignore TODO: @moritz error after I added the correct typing for DEFAULT_CIME_APP_CONFIG
        view: LineUpContext,
      },
    ],
  }
}

// TODO: @moritz We are currently using the global object directly, ideally we make it passable as prop.
export const CIMEAppContext = new API<AppState>(undefined, createRootReducer(CIMEReducers));

export function CIMEApp({ providePSEContext = true, ...props }: CIMEAppProps) {
  const [merged, setMerged] = React.useState<CIMEAppProps | null>(null);

  React.useEffect(() => {
    setMerged(merge(cloneDeep(DEFAULT_CIME_APP_CONFIG), props))
  }, [props.config, props.features, props.overrideComponents]);

  const component = <div style={{ width: '100%', height: '100%' }}>
    <Application
      ref={props.pseRef}
      config={merged?.config}
      features={merged?.features}
      //@ts-ignore
      overrideComponents={merged?.overrideComponents}
    />
  </div>


  return (providePSEContext ?
    <PSEContextProvider context={CIMEAppContext}>
      {component}
    </PSEContextProvider> : component
  );
}