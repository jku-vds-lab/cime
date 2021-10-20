import { useState } from "react";
import {
  PSEContextProvider,
  API,
  Application,
  PluginRegistry,
  createRootReducer,
} from "projection-space-explorer";
import { ChemPlugin } from "./Cime/ChemPlugin";
import { DatasetTabPanel } from "./Overrides/DatasetTabPanel";
import { CimeAppBar } from "./Overrides/CimeAppBar";
import { LineUpContext } from "./LineUpContext";
import { LineUpTabPanel } from "./Overrides/LineUpTabPanel";
import { AppState, CIMEReducers } from "./State/Store";

export const DEMO = false;

PluginRegistry.getInstance().registerPlugin(new ChemPlugin());

export function CIMEApp() {
  const [context] = useState(
    new API<AppState>(null, createRootReducer(CIMEReducers))
  );

  return (
    <PSEContextProvider context={context}>
      <Application
        config={{
          baseUrl: "http://localhost:5000",
          preselect: { url: "datasets/test.sdf" },
        }}
        features={{
          disableEmbeddings: {
            tsne: true,
            forceatlas: true,
          },
        }}
        overrideComponents={{
          datasetTab: DatasetTabPanel,
          appBar: CimeAppBar,
          tabs: [
            {
              name: "lineup",
              //@ts-ignore
              tab: LineUpTabPanel,
              title: "LineUp Integration",
              description: "Settings for LineUp Integration",
              icon: null,
            },
          ],
          detailViews: [
            {
              name: "lineup",
              //@ts-ignore
              view: LineUpContext,
            },
          ],
        }}
      />
    </PSEContextProvider>
  );
}
